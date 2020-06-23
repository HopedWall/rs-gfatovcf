use std::collections::HashMap;
use std::collections::BTreeMap; //like hashmap but sorted
use std::path::PathBuf;
use handlegraph::graph::HashGraph;
use handlegraph::graph::PathId;
use handlegraph::handlegraph::HandleGraph;
use handlegraph::handle::{NodeId,Direction,Handle};
use handlegraph::pathgraph::PathHandleGraph;
use std::fs::File;
use std::io::Write;
extern crate clap;
use clap::{Arg, App};
extern crate chrono;
use chrono::Utc;
use std::collections::VecDeque;
use std::collections::HashSet;
use std::cmp;
use std::cmp::Ordering;


struct Variant {
    chromosome: String,
    position: String,
    id: String,
    reference: String,
    alternate: String,
    quality: String,
    filter: String,
    info: String,
    format: String,
    sample_name: String
}


/// Returns a step as a String with NodeId and Orientation
fn process_step(h : &Handle) -> String {
    let is_rev = h.is_reverse();
    let id = h.id();
    
    //Added display trait in rs-handlegraph
    let mut result = id.to_string();
    
    if is_rev == true {
        result.push_str("+");
    } else {
        result.push_str("-");
    }

    result
}

/// Returns all paths as a hashmap, having the path_name as key and a list of steps as values
fn create_into_hashmap(g : &HashGraph, path_to_steps : &mut HashMap<String, Vec<String>>,  path : &PathId, step : &Handle) -> bool {    
    let path_name = g.get_path_name(path);
    if !path_to_steps.contains_key(path_name) {
        path_to_steps.insert(String::from(path_name), Vec::new());
    }
    path_to_steps.get_mut(path_name).unwrap().push(process_step(step));
    true
}

/// Prints an edge of a given HashGraph
fn show_edge(g_dfs : &HashGraph, a : &Handle, b : &Handle) {
    println!("{} --> {}", g_dfs.get_id(a), g_dfs.get_id(b));
}
/// Prints all nodes and edges of a given HashGraph
fn display_node_edges(g_dfs : &HashGraph, h : &Handle) {
   println!("node {}", g_dfs.get_id(h));
   g_dfs.follow_edges(
       h, Direction::Right,
       |n| {show_edge(&g_dfs,h,n); true}
   );
}

/// Calculates the distance of adjacent nodes
fn calculate_distance(visited_node_id_set : &mut HashSet<NodeId>, prev_node_id : &NodeId, neighbour_id : &NodeId, q : &mut VecDeque<NodeId>, distances_map : &mut BTreeMap<NodeId,u64>) {
    
    if !visited_node_id_set.contains(&neighbour_id) {
        //Clone necessary due to https://github.com/rust-lang/rust/issues/59159
        let previous_value = distances_map.get(prev_node_id).unwrap().clone();
        distances_map.insert(*neighbour_id, previous_value + 1);
        q.push_back(*neighbour_id);
        visited_node_id_set.insert(*neighbour_id);
    }
}

/// Finds the distance of each node from a given root
fn bfs_distances(g : &HashGraph, starting_node_id : &NodeId) -> (BTreeMap<NodeId,u64>, Vec<NodeId>) {
    let mut visited_node_id_set : HashSet<NodeId> = HashSet::new();
    let mut ordered_node_id_list : Vec<NodeId> = Vec::new();
    
    let mut distances_map : BTreeMap<NodeId, u64> = BTreeMap::new();
    let mut node_id_list = Vec::new();
    g.for_each_handle(|h| {
        node_id_list.push(g.get_id(h));
        true
    });
    for node_id in node_id_list {
        distances_map.insert(node_id, 0);
    }
    //node_id_list.clear();

    let mut q : VecDeque<NodeId> = VecDeque::new();
    q.push_back(*starting_node_id);
    visited_node_id_set.insert(*starting_node_id);
    while !&q.is_empty() {
        let current_node_id = q.pop_front().unwrap();
        let current_node = g.get_handle(current_node_id, false);
        ordered_node_id_list.push(current_node_id);

        g.follow_edges(
            &current_node, 
            Direction::Right, 
            |neighbor| {
                calculate_distance(&mut visited_node_id_set, 
                                   &current_node_id, &g.get_id(neighbor), 
                                   &mut q, 
                                   &mut distances_map);
                true
            });
    }

    (distances_map, ordered_node_id_list)
}

fn detect_bubbles(distances_map : &BTreeMap<NodeId,u64>, ordered_node_id_list : &Vec<NodeId>, 
                  dist_to_num_nodes : &BTreeMap<u64,usize>) -> Vec<(NodeId,NodeId)> {
    
    let mut possible_bubbles_list : Vec<(NodeId,NodeId)> = Vec::new();
    let mut first_bubble = true;

    for node_id in ordered_node_id_list {
        let mut pair : (NodeId,NodeId) = (NodeId::from(0),NodeId::from(0));
        let key = distances_map[&node_id];
        if dist_to_num_nodes[&key] == 1 {
            if !first_bubble {
                let latest_bubble = possible_bubbles_list.last_mut().unwrap();
                latest_bubble.1 = *node_id;
            }
            first_bubble = false;
            pair.0 = *node_id;
            possible_bubbles_list.push(pair); //Ok, pair.2 is a placeholder
        }
    }

    //Delete last bubble, won't be used anyway
    possible_bubbles_list.pop();

    possible_bubbles_list
}

//Computes a different bfs
fn new_bfs(g : &HashGraph, g_bfs : &mut HashGraph, node_id : &NodeId) {
    let current_handle = g.get_handle(*node_id, false);
    let mut added_handles : Vec<Handle> = vec![];

    if !g_bfs.has_node(*node_id) {
        g_bfs.create_handle(g.get_sequence(&current_handle), *node_id);
    }

    g.follow_edges(
        &current_handle,
        Direction::Right,
        |neighbor| {
            if !g_bfs.has_node(g.get_id(neighbor)) {
                let h = g_bfs.create_handle(
                    g.get_sequence(&neighbor), 
                    g.get_id(&neighbor)
                );
                added_handles.push(h);
                g_bfs.create_edge(&current_handle, neighbor);
            }
            
            true
        });
    
    for h in &added_handles {
        new_bfs(g, g_bfs, &g.get_id(h));
    }
        
}

/// Prints all paths between two nodes
fn print_all_paths_util(g : &HashGraph, u : &NodeId, d : &NodeId, visited_node_id_set : &mut HashSet<NodeId>, path_list : &mut Vec<NodeId>, all_path_list : &mut Vec<Vec<NodeId>>) {
    if !visited_node_id_set.contains(&u) {
        visited_node_id_set.insert(*u);
        path_list.push(*u);
    }

    //What does == mean exactly?
    if u == d {
        //TODO: check this
        all_path_list.push(path_list.to_vec());
    } else {
        g.follow_edges(
            &g.get_handle(*u, false), 
            Direction::Right,   //Is this correct?
            |i_node| {
                print_all_paths_util(g, &g.get_id(i_node), d, visited_node_id_set, path_list, all_path_list);
                true
            });
    }

    path_list.pop();
    visited_node_id_set.remove(&u);
}

/// Prints all paths in a given HashGrap, starting from a specific node and ending in another node
fn print_all_paths(g : &HashGraph, start_node_id : &NodeId, end_node_id : &NodeId, all_path_list : &mut Vec<Vec<NodeId>>) {
    let mut visited_node_id_set : HashSet<NodeId> = HashSet::new();
    
    let mut path_list : Vec<NodeId> = Vec::new();

    print_all_paths_util(g, start_node_id, end_node_id, &mut visited_node_id_set, &mut path_list, all_path_list);
}

fn detect_all_variants(path_to_steps_map : &HashMap<String,Vec<String>>, 
                       possible_bubbles_list : &Vec<(NodeId,NodeId)>,
                       graph : &HashGraph,
                       node_id_to_path_and_pos_map : &BTreeMap<NodeId, HashMap<String, usize>>,
                    ) -> Vec<Variant> {
    
    let mut stuff_to_alts_map : HashMap<String, HashSet<String>> = HashMap::new();
    
    // Taking each known path as reference, explore all bubbles in order to find variants;
    // these will be stored in stuff_to_alts_map
    for current_ref in path_to_steps_map.keys() {
        
        //Obtain all steps for each path
        let mut ref_path = Vec::new();
        println!("path_to_steps_map: {:?}",path_to_steps_map);
        for x in &path_to_steps_map[current_ref] {
            ref_path.push(x.parse::<u64>().unwrap());
        }

        detect_variants_per_reference(&current_ref,
                                      &ref_path, 
                                      possible_bubbles_list, 
                                      graph,
                                      node_id_to_path_and_pos_map,
                                      &mut stuff_to_alts_map);
    
    }

    // Convert stuff_to_alts_map to a more readable format
    let mut vcf_list : Vec<Variant> = Vec::new();
    for (chrom_pos_ref, alt_type_set) in &stuff_to_alts_map {
         
         let vec: Vec<&str> = chrom_pos_ref.split("_").collect();
         let chrom = vec[0];
         let pos = vec[1];
         let refr = vec[2];
         
         let mut alt_list : Vec<String> = Vec::new();
         for x in alt_type_set {
            let split : Vec<&str> = x.split("_").collect();
            alt_list.push(String::from(split[0]));
         }
         
         let mut type_set : Vec<&str> = Vec::new();
         for x in alt_type_set {
            let split : Vec<&str> = x.split("_").collect();
            type_set.push(split[1]);
         }

         let alts = alt_list.join(",");
         let types = type_set.join(";TYPE=");

         let v = Variant {
            chromosome: chrom.to_string(),
            position: pos.to_string(),
            id: ".".to_string(),
            reference: refr.to_string(),
            alternate: alts,
            quality: ".".to_string(),
            filter: ".".to_string(),
            info: format!("TYPE={}", types),
            format: "GT".to_string(),
            sample_name: "0|1".to_string(),            
        };

         vcf_list.push(v);
    }

    // Sort vcf_list for printing variants in the correct order
    vcf_list.sort_by(|a, b| {
        match a.chromosome.cmp(&b.chromosome) {
            Ordering::Equal => a.position.parse::<i32>().unwrap().cmp(&b.position.parse::<i32>().unwrap()),
            other => other,
        }});

    vcf_list
}
fn detect_variants_per_reference(current_ref : &String,
                                 ref_path : &Vec<u64>, 
                                 possible_bubbles_list : &Vec<(NodeId, NodeId)>,
                                 graph : &HashGraph,
                                 node_id_to_path_and_pos_map : &BTreeMap<NodeId, HashMap<String, usize>>,
                                 stuff_to_alts_map : &mut HashMap<String, HashSet<String>>
                                 ) {

    // Check all bubbles
    for (start,end) in possible_bubbles_list {
                
        println!("ref_path: {:?}",ref_path);
        println!("Bubble [{},{}]",start, end);

        let start_node_index_in_ref_path : usize;
        match ref_path.iter().position(|&r| NodeId::from(r) == *start) {
            None => continue,   //ignore, start not found in ref path
            Some(r) => start_node_index_in_ref_path = r,
        };

        let mut all_path_list : Vec<Vec<NodeId>> = Vec::new();
        print_all_paths(&graph, start, end, &mut all_path_list);

        //println!("All paths list: {:?}",all_path_list);
        for path in &all_path_list {
            println!("\tPath: {:?}", path);
            let mut pos_ref = node_id_to_path_and_pos_map[start][current_ref]+1;
            let mut pos_path = pos_ref;

            //println!("Start paths position: {}",pos_ref);

            let max_index = cmp::min(path.len(), ref_path.len());

            let mut current_index_step_path = 0;
            let mut current_index_step_ref = 0;

            for _i in 0..max_index {
                
                let mut current_node_id_ref = NodeId::from(ref_path[current_index_step_ref + start_node_index_in_ref_path]);
                let mut current_node_id_path = NodeId::from(path[current_index_step_path]);
                
                println!("{} {} ---> {} {}", pos_ref, pos_path, current_node_id_ref, current_node_id_path);
                if current_node_id_ref == current_node_id_path {
                    println!("REFERENCE");
                    let node_seq = graph.get_sequence(&graph.get_handle(current_node_id_ref, false));
                    pos_ref += node_seq.len();
                    pos_path = pos_ref;
                    
                    current_index_step_ref += 1;
                    current_index_step_path += 1;
                } else {

                    // Shouldn't be happening anymore
                    if current_index_step_path+1 >= path.len() {
                        break;
                    }
                    if current_index_step_ref + start_node_index_in_ref_path +1 >= ref_path.len() {
                        break;
                    }

                    let succ_node_id_path = NodeId::from(path[current_index_step_path + 1]);
                    let succ_node_id_ref = NodeId::from(ref_path[current_index_step_ref + start_node_index_in_ref_path + 1]);
                    if succ_node_id_ref == current_node_id_path {
                        
                        println!("DEL");
                        let node_seq_ref = graph.get_sequence(&graph.get_handle(current_node_id_ref, false));
                        
                        let prec_node_id_ref = NodeId::from(ref_path[current_index_step_ref + start_node_index_in_ref_path - 1]);
                        let prec_nod_seq_ref = graph.get_sequence(&graph.get_handle(prec_node_id_ref, false));
                        

                        let mut last = prec_nod_seq_ref.chars().last().unwrap().to_string();
                        //This must be concatenated!
                        last.push_str(&String::from(node_seq_ref));

                        let key = [current_ref.to_string(), (pos_path - 1).to_string(), last].join("_");
                        if !stuff_to_alts_map.contains_key(&key) {
                            stuff_to_alts_map.insert(key, HashSet::new());
                        }
                        //TODO: find a better way to do this
                        let mut last = prec_nod_seq_ref.chars().last().unwrap().to_string();
                        //This must be concatenated!
                        last.push_str(&String::from(node_seq_ref));
                        
                        let key = [current_ref.to_string(), (pos_path - 1).to_string(), last].join("_");
                        
                        let last = prec_nod_seq_ref.chars().last().unwrap().to_string();
                        let mut string_to_insert = last;
                        string_to_insert.push_str("_del");
                        stuff_to_alts_map.get_mut(&key).unwrap().insert(string_to_insert);
                        
                        pos_ref += node_seq_ref.len();

                        current_index_step_ref += 1;
                        current_node_id_ref = NodeId::from(ref_path[current_index_step_ref + start_node_index_in_ref_path -1]);
                        println!("\t {}", current_node_id_ref);
                        continue;
                    } else if succ_node_id_path == current_node_id_ref {
                        println!("INS");
                        let node_seq_path = graph.get_sequence(&graph.get_handle(current_node_id_path,false));
                        
                        let prec_node_id_ref = NodeId::from(ref_path[current_index_step_ref + start_node_index_in_ref_path-1]);
                        let prec_nod_seq_ref = graph.get_sequence(&graph.get_handle(prec_node_id_ref,false));

                        let last = prec_nod_seq_ref.chars().last().unwrap().to_string();
                        //let key = [current_ref.to_string(), (pos_ref-1).to_string(), String::from(prec_nod_seq_ref)].join("_");
                        let key = [current_ref.to_string(), (pos_ref-1).to_string(), last].join("_");
                        if !stuff_to_alts_map.contains_key(&key) {
                            stuff_to_alts_map.insert(key, HashSet::new());
                        }
                        
                        //Re-create key since it goes out of scope
                        let last = prec_nod_seq_ref.chars().last().unwrap().to_string();
                        let key = [current_ref.to_string(), (pos_ref-1).to_string(), last].join("_");
                        
                        let last = prec_nod_seq_ref.chars().last().unwrap().to_string();
                        let mut string_to_insert = String::from(last);
                        string_to_insert.push_str(node_seq_path);
                        string_to_insert.push_str("_ins");
                        stuff_to_alts_map.get_mut(&key).unwrap().insert(string_to_insert);

                        pos_path += node_seq_path.len();
                        
                        current_index_step_path += 1;
                        current_node_id_path = path[current_index_step_path];
                        println!("\t{}", current_node_id_path);
                        continue;
                    } else {
                        let node_seq_ref = graph.get_sequence(&graph.get_handle(current_node_id_ref, false));
                        let node_seq_path = graph.get_sequence(&graph.get_handle(current_node_id_path, false));

                        if node_seq_ref == node_seq_path {
                            println!("REFERENCE");
                        } else {
                            println!("SNV");
                        }

                        let key = [current_ref.to_string(), pos_path.to_string(), node_seq_ref.to_string()].join("_");
                        if !stuff_to_alts_map.contains_key(&key){
                            stuff_to_alts_map.insert(key, HashSet::new());
                        }
                                                        
                        //TODO: find a better way to do this
                        let key = [current_ref.to_string(), pos_path.to_string(), node_seq_ref.to_string()].join("_");
                        let mut string_to_insert = String::from(node_seq_path.chars().last().unwrap().to_string());
                        string_to_insert.push_str("_snv");
                        stuff_to_alts_map.get_mut(&key).unwrap().insert(string_to_insert);


                        pos_ref += node_seq_ref.len();
                        pos_path += node_seq_path.len();
                        current_index_step_ref += 1;
                        current_index_step_path += 1;
                    }
                }
            }
            println!("---");
        }   
    }
    println!("==========================================");
}

fn get_node_positions_in_paths(graph : &HashGraph, path_to_steps_map : &mut HashMap<String,Vec<String>>) -> BTreeMap<NodeId, HashMap<String, usize>> {
    let mut node_id_to_path_and_pos_map : BTreeMap<NodeId, HashMap<String, usize>> = BTreeMap::new();

    for (path_name, steps_list) in path_to_steps_map {
        let mut pos = 0;

        for node_id_is_rev in steps_list {
            
            // Get orientation
            let _is_rev = node_id_is_rev.pop().unwrap();
            // Get the id of the node string -> NodeId
            let node_id : NodeId = NodeId::from(node_id_is_rev.parse::<u64>().unwrap());
            
            let node_handle = graph.get_handle(node_id, false);
            let seq = graph.get_sequence(&node_handle);

            if !node_id_to_path_and_pos_map.contains_key(&node_id) {
                node_id_to_path_and_pos_map.insert(node_id, HashMap::new());
            }

            if !node_id_to_path_and_pos_map[&node_id].contains_key(path_name) {
                node_id_to_path_and_pos_map.get_mut(&node_id).unwrap().insert(String::from(path_name), pos);
            }

            pos += seq.len();
        }
    }
    
    node_id_to_path_and_pos_map
}

fn paths_to_steps(graph : &HashGraph) -> HashMap<String,Vec<String>> {
    let mut path_to_steps_map : HashMap<String,Vec<String>> = HashMap::new();

    graph.for_each_path_handle(|p| {
        graph.for_each_step_in_path(&p, 
                                    |s| {
                                        create_into_hashmap(&graph, 
                                                            &mut path_to_steps_map, 
                                                            &p, 
                                                            &s);
                                        true
                                        });
        true
    });

    path_to_steps_map
}

fn get_dist_to_num_nodes(distances_map : &BTreeMap<NodeId, u64>) -> BTreeMap<u64,usize> {
    
    let mut dist_to_num_nodes : BTreeMap<u64,usize> = BTreeMap::new();

    for (_, distance) in distances_map.iter() {
        if !dist_to_num_nodes.contains_key(&distance) {
            dist_to_num_nodes.insert(*distance,0);
        }
        *dist_to_num_nodes.get_mut(distance).unwrap() += 1;
    }

    dist_to_num_nodes
}

fn get_path_to_sequence(graph : &HashGraph, path_to_steps_map : &HashMap<String,Vec<String>>) -> HashMap<String,String> {
    let mut path_to_sequence_map : HashMap<String,String> = HashMap::new();

    for (path_name, steps_list) in path_to_steps_map {
        path_to_sequence_map.insert(path_name.to_string(), String::new());

         for node_id_rev in steps_list {
             path_to_sequence_map.get_mut(path_name).unwrap().push_str(graph.get_sequence(&graph.get_handle(NodeId::from(node_id_rev.parse::<u64>().unwrap()), false)));
         }
    }

    path_to_sequence_map
}

fn main() {

    let matches = App::new("rs-GFAtoVCF")
                          .version("1.0")
                          .author("Francesco Porto <francesco.porto97@gmail.com>")
                          .about("Converts GFA to VCF")
                          .arg(Arg::with_name("input")
                               .short("i")
                               .long("input")
                               .value_name("FILE")
                               .help("Sets the input file to use")
                               .required(true)
                               .takes_value(true))
                          .arg(Arg::with_name("output")
                               .short("o")
                               .long("output")
                               .value_name("FILE")
                               .help("Sets the output file to use")
                               .required(true)
                               .takes_value(true))
                          .get_matches();

    let in_path_file = matches.value_of("input").expect("Could not parse argument --input");
    let out_path_file = matches.value_of("output").expect("Could not parse argument --output");

    //let in_path_file = "./input/samplePath3.gfa";
    //let out_path_file = "./input/samplePath3.vcf";
      
    if let Some(gfa) = gfa::parser::parse_gfa(&PathBuf::from(in_path_file)) {
        
        let graph = HashGraph::from_gfa(&gfa);

        // Obtains, for each path, a list of all its steps (its nodes)
        let mut path_to_steps_map : HashMap<String,Vec<String>> = paths_to_steps(&graph);
        
        // Obtains, for each node, its position in each path where the node is in
        let node_id_to_path_and_pos_map : BTreeMap<NodeId, HashMap<String, usize>> = get_node_positions_in_paths(&graph, &mut path_to_steps_map);

        // for node_id in node_id_to_path_and_pos_map.keys() {
        //     let path_and_pos_map = node_id_to_path_and_pos_map.get(node_id); 
        //     println!("Node_id : {}", node_id);

        //     for (path, pos) in path_and_pos_map.unwrap() {
        //         println!("Path: {}  -- Pos: {}", path, pos);
        //     }
        // }

        let mut g_bfs = HashGraph::new();

        new_bfs(&graph, &mut g_bfs, &NodeId::from(1));

        // g_bfs.for_each_handle(|h| {
        //                             display_node_edges(&g_bfs, &h);
        //                             true
        //                           });
        

        // Obtains, for each level of the tree, how many nodes are there
        let (distances_map, ordered_node_id_list) = bfs_distances(&g_bfs, &NodeId::from(1)); 
        
        // println!("\nNode --> Distance from root");
        // for (node_id, distance) in distances_map.iter() {
        //     println!("{} - distance from root: {}", node_id, distance);
        // }

        //Obtain a map where, for each distance from root, the number of nodes
        //at that distance are present
        let dist_to_num_nodes : BTreeMap<u64,usize> = get_dist_to_num_nodes(&distances_map);
        
        println!("\nDistance from root --> Num. nodes");
        for (k,v) in dist_to_num_nodes.iter() {
            println!("{} --> {}", k, v);
        }

        println!("\nBubbles");
        let possible_bubbles_list : Vec<(NodeId,NodeId)> = detect_bubbles(&distances_map, &ordered_node_id_list, &dist_to_num_nodes);

        println!("Detected bubbles {:#?}",possible_bubbles_list);
        println!("\n------------------");

        //Obtains, for each path, a string representing all the bases of the path (not actually used)
        //let _path_to_sequence_map : HashMap<String,String> = get_path_to_sequence(&graph, &path_to_steps_map);
        //println!("Path to sequence: {:?}",path_to_sequence_map);

        let vcf_list = detect_all_variants(&path_to_steps_map, 
                                           &possible_bubbles_list, 
                                           &graph, 
                                           &node_id_to_path_and_pos_map);

        //Write variants to file
        write_to_file(&PathBuf::from(out_path_file), &vcf_list).unwrap();
    } else {
        panic!("Couldn't parse gfa file!");
    }

}

fn write_to_file(path: &PathBuf, variations: &Vec<Variant>) -> std::io::Result<()> {
    let mut file = File::create(path).expect(&format!("Error creating file {:?}", path));

    let header = [
        "##fileformat=VCFv4.2",
        &format!("##fileDate={}",Utc::now().format("%Y-%m-%d %H:%M:%S").to_string()),
        "##reference=x.fa",
        "##reference=y.fa",
        "##reference=z.fa",
        "##INFO=<ID=TYPE,Number=A,Type=String,Description=\"Type of each allele (snv, ins, del, mnp, complex)\">",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        &["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SampleName"].join("\t"),
    ].join("\n");
    file.write(header.as_bytes()).expect("Error writing header");

    file.write("\n".as_bytes()).expect("Error writing to file");
    for var in variations {
        let to_write = format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                                var.chromosome, 
                                var.position, 
                                var.id, 
                                var.reference, 
                                var.alternate, 
                                var.quality, 
                                var.filter,
                                var.info,
                                var.format,
                                var.sample_name);
        file.write(to_write.as_bytes()).expect("Error writing variant");
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use handlegraph::graph::HashGraph;
    use handlegraph::handlegraph::HandleGraph;
    use super::*;

    //Used in other tests
    fn read_test_gfa() -> HashGraph {
        use gfa::parser::parse_gfa;

        HashGraph::from_gfa(&parse_gfa(&PathBuf::from("./input/samplePath3.gfa")).unwrap())
    }

    #[test]
    fn test_path_to_steps() {
        
        let graph = read_test_gfa();

        let mut path_to_steps_map : HashMap<String,Vec<String>> = HashMap::new();

        graph.for_each_path_handle(|p| {
            graph.for_each_step_in_path(&p, 
                                        |s| {
                                            create_into_hashmap(&graph, 
                                                                &mut path_to_steps_map, 
                                                                &p, 
                                                                &s);
                                            true
                                            });
            true
        });

        // Check if all paths have been found
        assert_eq!(path_to_steps_map.keys().len(), 3);
        assert!(path_to_steps_map.contains_key("x"));
        assert!(path_to_steps_map.contains_key("y"));
        assert!(path_to_steps_map.contains_key("z"));
        

        //Check for each path that all its node have been found

        //P	x	1+,3+,5+,6+,8+,9+,11+,12+,13+,15+,16+,18+,19+
        let path_x : Vec<String> = vec!["1+".to_string(),
                                        "3+".to_string(),
                                        "5+".to_string(),
                                        "6+".to_string(),
                                        "8+".to_string(),
                                        "9+".to_string(),
                                        "11+".to_string(),
                                        "12+".to_string(),
                                        "13+".to_string(),
                                        "15+".to_string(),
                                        "16+".to_string(),
                                        "18+".to_string(),
                                        "19+".to_string()];
        assert_eq!(*path_to_steps_map.get("x").unwrap(), path_x);

        //P	y	1+,2+,5+,6+,8+,9+,10+,11+,13+,14+,16+,17+,19+
        let path_y : Vec<String> = vec!["1+".to_string(),
                                        "2+".to_string(),
                                        "5+".to_string(),
                                        "6+".to_string(),
                                        "8+".to_string(),
                                        "9+".to_string(),
                                        "10+".to_string(),
                                        "11+".to_string(),
                                        "13+".to_string(),
                                        "14+".to_string(),
                                        "16+".to_string(),
                                        "17+".to_string(),
                                        "19+".to_string()];
        assert_eq!(*path_to_steps_map.get("y").unwrap(), path_y);

        //P	z	1+,2+,4+,6+,7+,9+,10+,11+,12+,13+,14+,16+,18+,19+
        let path_z : Vec<String> = vec!["1+".to_string(),
                                        "2+".to_string(),
                                        "4+".to_string(),
                                        "6+".to_string(),
                                        "7+".to_string(),
                                        "9+".to_string(),
                                        "10+".to_string(),
                                        "11+".to_string(),
                                        "12+".to_string(),
                                        "13+".to_string(),
                                        "14+".to_string(),
                                        "16+".to_string(),
                                        "18+".to_string(),
                                        "19+".to_string()];
        assert_eq!(*path_to_steps_map.get("z").unwrap(), path_z);        
    }

    #[test]
    fn test_bfs() {
        let graph = read_test_gfa();
        let mut g_bfs = HashGraph::new();
        new_bfs(&graph, &mut g_bfs, &NodeId::from(1));

        // All nodes must be present in bfs
        assert_eq!(graph.get_node_count(), g_bfs.get_node_count());

        // There should be less (or the same number of) edges in g_bfs
        // since nodes can only get added once
        assert!(graph.get_edge_count() >= g_bfs.get_edge_count());
    }

    #[test]
    fn test_bfs_distances() {
        let graph = read_test_gfa();
        let mut g_bfs = HashGraph::new();
        new_bfs(&graph, &mut g_bfs, &NodeId::from(1));
        
        let (distances_map, ordered_node_id_list) = bfs_distances(&g_bfs, &NodeId::from(1));

        // Should pass easily
        assert_eq!(distances_map[&NodeId::from(2)], distances_map[&NodeId::from(3)]);
        assert_eq!(distances_map[&NodeId::from(4)], distances_map[&NodeId::from(5)]);
        assert_eq!(distances_map[&NodeId::from(7)], distances_map[&NodeId::from(8)]);

        // Caused problems in dfs
        assert_eq!(distances_map[&NodeId::from(10)], distances_map[&NodeId::from(11)]);
        assert_eq!(distances_map[&NodeId::from(12)], distances_map[&NodeId::from(13)]);

        //Should pass easily
        assert_eq!(distances_map[&NodeId::from(14)], distances_map[&NodeId::from(15)]);
        assert_eq!(distances_map[&NodeId::from(17)], distances_map[&NodeId::from(18)]);

        // Check ordered_node_list now
        let mut sorted = ordered_node_id_list.clone();
        sorted.sort();
        assert!(ordered_node_id_list.eq(&sorted));
    }

    #[test]
    fn test_dist_to_num_nodes() {
        let graph = read_test_gfa();
        let mut g_bfs = HashGraph::new();
        new_bfs(&graph, &mut g_bfs, &NodeId::from(1));
        
        let (distances_map, _) = bfs_distances(&g_bfs, &NodeId::from(1));
        let mut dist_to_num_nodes : BTreeMap<u64,usize> = BTreeMap::new();
        for (_, distance) in distances_map.iter() {
            if !dist_to_num_nodes.contains_key(&distance) {
                dist_to_num_nodes.insert(*distance,0);
            }
            *dist_to_num_nodes.get_mut(distance).unwrap() += 1;
        }

        //println!("{:#?}",dist_to_num_nodes);

        //Root
        assert_eq!(dist_to_num_nodes[&0], 1);

        assert_eq!(dist_to_num_nodes[&1], 2);
        assert_eq!(dist_to_num_nodes[&2], 2);
        assert_eq!(dist_to_num_nodes[&3], 1);
        assert_eq!(dist_to_num_nodes[&4], 2);
        assert_eq!(dist_to_num_nodes[&5], 1);

        // Critical, did not work with dfs
        assert_eq!(dist_to_num_nodes[&6], 2);
        assert_eq!(dist_to_num_nodes[&7], 2);
        assert_eq!(dist_to_num_nodes[&8], 2);

        assert_eq!(dist_to_num_nodes[&9], 1);
        assert_eq!(dist_to_num_nodes[&10], 2);
        assert_eq!(dist_to_num_nodes[&11], 1);
    }

    #[test]
    fn test_bubble_detection() {
        let graph = read_test_gfa();
        let mut g_bfs = HashGraph::new();
        new_bfs(&graph, &mut g_bfs, &NodeId::from(1));
        
        let (distances_map, ordered_node_id_list) = bfs_distances(&g_bfs, &NodeId::from(1));
        let mut dist_to_num_nodes : BTreeMap<u64,usize> = BTreeMap::new();
        for (_, distance) in distances_map.iter() {
            if !dist_to_num_nodes.contains_key(&distance) {
                dist_to_num_nodes.insert(*distance,0);
            }
            *dist_to_num_nodes.get_mut(distance).unwrap() += 1;
        }

        let possible_bubbles_list : Vec<(NodeId,NodeId)> = detect_bubbles(&distances_map, &ordered_node_id_list, &dist_to_num_nodes);

        //println!("Possible bubbles list: {:#?}",possible_bubbles_list);
        
        assert!(possible_bubbles_list.contains(&(NodeId::from(1), NodeId::from(6))));
        assert!(possible_bubbles_list.contains(&(NodeId::from(6), NodeId::from(9))));
        assert!(possible_bubbles_list.contains(&(NodeId::from(9), NodeId::from(16))));
        assert!(possible_bubbles_list.contains(&(NodeId::from(16), NodeId::from(19))));
    }

    fn run_whole_script() {
        
    }

    #[test]
    fn test_variant_detection() {
        //TODO: think how to easily do it!
    }

}