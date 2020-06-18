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
//use clap::{Arg, App};
extern crate chrono;
use chrono::Utc;
use std::io::{BufReader,BufRead};
use std::collections::VecDeque;
use std::collections::HashSet;
use std::cmp;
use std::iter::FromIterator;
use std::cmp::Ordering;


struct Variant {
    chromosome: String,
    position: i32,
    id: String,
    reference: String,
    alternate: String,
    quality: i32,
    filter: String,
    info: String,
}


/// Returns a step as a String with NodeId and Orientation
fn process_step(g : &HashGraph, h : &Handle) -> String {
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
    path_to_steps.get_mut(path_name).unwrap().push(process_step(g, step));
    true
}

/// Adds edges to DFS
fn create_edge_and_so_on(g : &HashGraph, g_dfs : &mut HashGraph, handle1 : &Handle, handle2 : &Handle, neighbor_node_id : &NodeId) {
    let handle1_id = g.get_id(handle1);
    let handle2_id = g.get_id(handle2); 
    //println!("Create edge from {} to {}",handle1_id,handle2_id);

    if !g_dfs.has_node(handle2_id) {
        //println!("g_dfs does not have node: {}",handle2_id);
        dfs(g,g_dfs,neighbor_node_id);
        
        if !g_dfs.has_node(handle1_id) {
            //println!("Add node {}",handle1_id);
            g_dfs.create_handle(
                g.get_sequence(handle1), 
                handle1_id);
        }
    
        if !g_dfs.has_node(handle2_id) {
            //println!("Add node {}",handle2_id);
            g_dfs.create_handle(
                g.get_sequence(handle2), 
                handle2_id);
        }
    
        g_dfs.create_edge(handle1, handle2)
    }
}
/// Computes the DFS of a given HashGraph
fn dfs(g : &HashGraph, g_dfs : &mut HashGraph, node_id : &NodeId) {
    let current_node = g.get_handle(*node_id, false);
    //let sequence_node = g.get_sequence(&current_node);

    g.follow_edges(
        &current_node, 
        Direction::Right,  //What should go here?
        |neighbor| {
            create_edge_and_so_on(
                &g, g_dfs, 
                &current_node, neighbor,
                &g.get_id(neighbor));
                true
        });
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
fn calculate_distance(visited_node_id_set : &mut HashSet<NodeId>, prev_node_id : &NodeId, neighbour_id : &NodeId, q : &mut VecDeque<NodeId>, distances_map : &mut HashMap<NodeId,u64>) {
    
    if !visited_node_id_set.contains(&neighbour_id) {
        let previous_value = distances_map.get(prev_node_id).unwrap();
        distances_map.insert(*neighbour_id, *previous_value + 1);
        q.push_back(*neighbour_id);
        visited_node_id_set.insert(*neighbour_id);
    }
}

// This function can be optimized via closures!
/// Finds the distance of each node from root
fn bfs_distances(g : &HashGraph, starting_node_id : &NodeId) -> (HashMap<NodeId,u64>, Vec<NodeId>) {
    let mut visited_node_id_set : HashSet<NodeId> = HashSet::new();
    let mut ordered_node_id_list : Vec<NodeId> = Vec::new();
    
    let mut distances_map : HashMap<NodeId, u64> = HashMap::new();
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
                calculate_distance(&mut visited_node_id_set, &current_node_id, &g.get_id(neighbor), &mut q, &mut distances_map);
                true
            });
    }

    (distances_map, ordered_node_id_list)
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

fn main() {

    // let matches = App::new("rs-GFAtoVCF")
    //                       .version("1.0")
    //                       .author("Francesco Porto <francesco.porto97@gmail.com>")
    //                       .about("Converts GFA to VCF")
    //                       .arg(Arg::with_name("input")
    //                            .short("i")
    //                            .long("input")
    //                            .value_name("FILE")
    //                            .help("Sets the input file to use")
    //                            .required(true)
    //                            .takes_value(true))
    //                       .arg(Arg::with_name("output")
    //                            .short("o")
    //                            .long("output")
    //                            .value_name("FILE")
    //                            .help("Sets the output file to use")
    //                            .required(true)
    //                            .takes_value(true))
    //                       .get_matches();

    // let in_path = matches.value_of("input").expect("Could not parse argument --input");
    // let out_path = matches.value_of("output").expect("Could not parse argument --output");

    let in_path_file = "./input/samplePath3.gfa";
    let out_path_file = "./input/samplePath3.vcf";
      
    if let Some(gfa) = gfa::parser::parse_gfa(&PathBuf::from(in_path_file)) {
        
        let graph = HashGraph::from_gfa(&gfa);
       // println!("{:?}",graph);

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

        println!("Path to steps map");
        println!("{:?}",path_to_steps_map);
        
        // Obtains, for each node_id, the starting position in terms of actual sequence length, according to path
        let mut node_id_to_path_and_pos_map : HashMap<NodeId, HashMap<String, usize>> = HashMap::new();
        for (path_name, steps_list) in &mut path_to_steps_map {
            let mut pos = 0;

            for node_id_is_rev in steps_list {
                
                // Get orientation
                let is_rev = node_id_is_rev.pop().unwrap();
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

        //println!("Node id to path and pos");
        //println!("{:?}",node_id_to_path_and_pos_map);

        // node_id_to_path_and_pos_map must be sorted
        // TODO: find a better inplace solution
        // maybe use btreemap directly
        let mut sorted = BTreeMap::new();
        for (id, value) in node_id_to_path_and_pos_map.iter() {
            sorted.insert(id, value);
        }

        //println!("Sorted");
        //println!("{:?}",sorted);

        for node_id in sorted.keys() {
            let path_and_pos_map = node_id_to_path_and_pos_map.get(node_id); 
            println!("Node_id : {}", node_id);

            for (path, pos) in path_and_pos_map.unwrap() {
                println!("Path: {}  -- Pos: {}", path, pos);
            }
        }

        let start_node = graph.get_handle(NodeId::from(1), false);

        let mut g_dfs = HashGraph::new();

        // Compute dfs
        dfs(&graph, &mut g_dfs, &NodeId::from(1));

        // let h1=g_dfs.create_handle("CAAATAAG", NodeId::from(1));
        // let h2=g_dfs.create_handle("A", NodeId::from(2));
        // let h3=g_dfs.create_handle("G", NodeId::from(3));
        // let h4=g_dfs.create_handle("T", NodeId::from(4));
        // let h5=g_dfs.create_handle("C", NodeId::from(5));
        // let h6=g_dfs.create_handle("TTG", NodeId::from(6));
        // let h7=g_dfs.create_handle("A", NodeId::from(7));
        // let h8=g_dfs.create_handle("G", NodeId::from(8));
        // let h9=g_dfs.create_handle("AAAT", NodeId::from(9));
        // let h10=g_dfs.create_handle("AA",NodeId::from(10));
        // let h11=g_dfs.create_handle("TTTCT", NodeId::from(11));
        // let h12=g_dfs.create_handle("GG",NodeId::from(12));
        // let h13=g_dfs.create_handle("AGTTCTAT", NodeId::from(13));
        // let h14=g_dfs.create_handle("A", NodeId::from(14));
        // let h15=g_dfs.create_handle("T", NodeId::from(15));
        // let h16=g_dfs.create_handle("ATAT", NodeId::from(16));
        // let h17=g_dfs.create_handle("A", NodeId::from(17));
        // let h18=g_dfs.create_handle("T", NodeId::from(18));
        // let h19=g_dfs.create_handle("CCAACTCTCTG", NodeId::from(19));

        // g_dfs.create_edge(&h1, &h2);
        // g_dfs.create_edge(&h1, &h3);
        // g_dfs.create_edge(&h3, &h4);
        // g_dfs.create_edge(&h3, &h5);
        // g_dfs.create_edge(&h5, &h6);
        // g_dfs.create_edge(&h6, &h7);
        // g_dfs.create_edge(&h6, &h8);
        // g_dfs.create_edge(&h8, &h9);
        // g_dfs.create_edge(&h9, &h10);
        // g_dfs.create_edge(&h9, &h11);
        // g_dfs.create_edge(&h11, &h12);
        // g_dfs.create_edge(&h11, &h13);
        // g_dfs.create_edge(&h13, &h14);
        // g_dfs.create_edge(&h13, &h15);
        // g_dfs.create_edge(&h15, &h16);
        // g_dfs.create_edge(&h16, &h17);
        // g_dfs.create_edge(&h16, &h18);
        // g_dfs.create_edge(&h18, &h19);

        // println!("g_dfs is: {:?}",g_dfs);
   

        g_dfs.for_each_handle(|h| {
                                    display_node_edges(&g_dfs, &h);
                                    true
                                  });
        
        let (distances_map, ordered_node_id_list) = bfs_distances(&g_dfs, &NodeId::from(1)); 
        //println!("distances map: {:?}",distances_map);
        //println!("ordered_node_id_list: {:?}",ordered_node_id_list);
        
        // TODO: find a better inplace solution
        // maybe use btreemap directly
        let mut sorted_2 = BTreeMap::new();
        for (id, value) in distances_map.iter() {
            sorted_2.insert(id, value);
        }
        println!("\nNode --> Distance from root");
        for (node_id, distance) in sorted_2.iter() {
            println!("{} - distance from root: {}", node_id, distance);
        }

        let mut dist_to_num_nodes : HashMap<u64,usize> = HashMap::new();

        for (_, distance) in distances_map.iter() {
            if !dist_to_num_nodes.contains_key(&distance) {
                dist_to_num_nodes.insert(*distance,0);
            }
            *dist_to_num_nodes.get_mut(distance).unwrap() += 1;
        }

        
        // TODO: find a better inplace solution
        // maybe use btreemap directly
        let mut sorted_3 = BTreeMap::new();
        for (id, value) in dist_to_num_nodes.iter() {
            sorted_3.insert(id, value);
        }
        
        println!("\nDistance from root --> Num. nodes");
        for (k,v) in sorted_3.iter() {
            println!("{} --> {}", k, v);
        }

        println!("\nBubbles");
        let mut possible_bubbles_list : Vec<(NodeId,NodeId)> = Vec::new();
        let mut first_bubble = true;

        for node_id in ordered_node_id_list {
            let mut pair : (NodeId,NodeId) = (NodeId::from(0),NodeId::from(0));
            let key = distances_map[&node_id];
            if dist_to_num_nodes[&key] == 1 {
                if !first_bubble {
                    println!("{} END {:?} {}", node_id, node_id_to_path_and_pos_map[&node_id],g_dfs.get_sequence(&g_dfs.get_handle(node_id, false)));
                    
                    let latest_bubble = possible_bubbles_list.last_mut().unwrap();
                    latest_bubble.1 = node_id;
                }
                first_bubble = false;
                println!("{} START {:?} {}",node_id, node_id_to_path_and_pos_map[&node_id],g_dfs.get_sequence(&g_dfs.get_handle(node_id, false)));
                pair.0 = node_id;
                possible_bubbles_list.push(pair); //Ok, pair.2 is a placeholder
                //println!("Possible bubbles: {:?}",possible_bubbles_list);
            } else {
                println!("{} Bubble {:?} {:?}",node_id, node_id_to_path_and_pos_map[&node_id], g_dfs.get_sequence(&g_dfs.get_handle(node_id,false)));
            }
        }

        println!("Possible bubbles list: {:?}",possible_bubbles_list);
        possible_bubbles_list.pop(); //Remove last bubble, should not be considered
        println!("Possible bubbles list: {:?}",possible_bubbles_list);

        //Compress bubble list
        let mut possible_bubbles_list_compressed : Vec<(NodeId,NodeId)> = Vec::new();
        let mut temp : (NodeId, NodeId) = (NodeId::from(0),NodeId::from(0));
        let mut open_bubble = false;
        for (start, end) in possible_bubbles_list {
            if end==start+1 {
                if !open_bubble {
                    temp.0 = start;
                    temp.1 = end;
                    open_bubble = true;
                } else {
                    temp.1 = end;
                }
            } else {
                if open_bubble {
                    
                    if temp.1 == start {
                        // Extend open bubble
                        temp.1 = end;

                        possible_bubbles_list_compressed.push(temp);
                        temp.0 = NodeId::from(0);
                        temp.1 = NodeId::from(0);
                        open_bubble = false;

                        continue;

                        //Do not push (start,end) since already inserted
                    }

                    possible_bubbles_list_compressed.push(temp);
                    temp.0 = NodeId::from(0);
                    temp.1 = NodeId::from(0);
                    open_bubble = false;
                }
                possible_bubbles_list_compressed.push((start, end));
            }
        }

        println!("Possible bubbles list compressed: {:?}",possible_bubbles_list_compressed);

        println!("\n------------------");

        let mut path_to_sequence_map : HashMap<String,String> = HashMap::new();
        for (path_name, steps_list) in &path_to_steps_map {
            path_to_sequence_map.insert(path_name.to_string(), String::new());

             for node_id_rev in steps_list {
                 path_to_sequence_map.get_mut(path_name).unwrap().push_str(graph.get_sequence(&graph.get_handle(NodeId::from(node_id_rev.parse::<u64>().unwrap()), false)));
             }
        }

        println!("Path to sequence: {:?}",path_to_sequence_map);

        let mut stuff_to_alts_map : HashMap<String, HashSet<String>>  = HashMap::new();
        
        //Consider each path as reference
        for current_ref in path_to_steps_map.keys() {
            //println!("LOOP START {:?}",current_ref);
            
            //Obtain all steps for each path
            //TODO: maybe rewrite in a more compact way
            let mut ref_path = Vec::new();
            println!("path_to_steps_map: {:?}",path_to_steps_map);
            for x in &path_to_steps_map[current_ref] {
                ref_path.push(x.parse::<u64>().unwrap());
            }
 
            // Evaluate each possible bubble on reference path
            for (start,end) in &possible_bubbles_list_compressed {
                
                println!("ref_path: {:?}",ref_path);
                println!("Bubble [{},{}]",start, end);
                //println!("Possible bubbles list {:?}",possible_bubbles_list);

                //let start_node_index_in_ref_path = ref_path.iter().position(|&r| NodeId::from(r) == *start).unwrap();

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

                    println!("Start paths position: {}",pos_ref);

                    //let max_index = if path.len() < ref_path.len() {path.len()} else {ref_path.len()};
                    let max_index = cmp::min(path.len(), ref_path.len());

                    //println!("\nDEBUG:");
                    //println!("Path: {:?}", path);
                    //println!("Ref_Path: {:?}\n", ref_path);
                    //println!("path: {} ref_path: {}", path.len(), ref_path.len());
                    //println!("max_index: {}",max_index);

                    let mut current_index_step_path = 0;
                    let mut current_index_step_ref = 0;

                    for _i in 0..max_index {
                        //println!("\nDEBUG:");
                        println!("Current index step path: {}",current_index_step_path);
                        println!("Current index step ref: {}",current_index_step_ref);
                        println!("Max index: {}",max_index);
                        println!("Start node index in ref path: {}\n",start_node_index_in_ref_path);
                        
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

                            if current_index_step_path+1 >= path.len() {
                                //TODO: figure out what this means
                                break;
                            }

                            if current_index_step_ref + start_node_index_in_ref_path +1 >= ref_path.len() {
                                //TODO: figure out what this means
                                break;
                            }

                            //Index out of bounds happens here
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

        println!("Stuff to alts map: {:#?}",stuff_to_alts_map);
        let mut vcf_list : Vec<Vec<String>> = Vec::new();
        for (chrom_pos_ref, alt_type_set) in &stuff_to_alts_map {
             
             let vec: Vec<&str> = chrom_pos_ref.split("_").collect();
             //print!("Vec: {:?}",vec);
             let chrom = vec[0];
             let pos = vec[1];
             let refr = vec[2];
             
             let mut alt_list : Vec<String> = Vec::new();
             for x in alt_type_set {
                //println!("x is {}",x);
                let split : Vec<&str> = x.split("_").collect();
                alt_list.push(String::from(split[0]));
             }
             
             //let type_set : HashSet<String> = HashSet::new(); 
             let mut type_set : Vec<&str> = Vec::new();
             for x in alt_type_set {
                let split : Vec<&str> = x.split("_").collect();
                //type_set.insert(String::from(split[1]));
                type_set.push(split[1]);
             }

             let alts = alt_list.join(",");
             let types = type_set.join(";TYPE=");
             let list_to_append = [chrom, pos,".",refr,&alts,".",".", "TYPE=",&types, "GT", "0|1"];
             vcf_list.push(Vec::from_iter(list_to_append.iter().map(|&x| String::from(x))));
        }

        //println!("VCF list is {:#?}",vcf_list);

        //sort by chrom
        //vcf_list.sort_by(|a, b| a.partial_cmp(b).unwrap());
        
        // First sort by chrom, then by starting pos
        vcf_list.sort_by(|a, b| {
            match a[0].cmp(&b[0]) {
                Ordering::Equal => a[1].parse::<i32>().unwrap().cmp(&b[1].parse::<i32>().unwrap()),
                other => other,
            }});


        //Write variants to file
        write_to_file(&PathBuf::from(out_path_file), &vcf_list).unwrap();
    } else {
        panic!("Couldn't parse gfa file!");
    }

}

fn write_to_file(path: &PathBuf, variations: &Vec<Vec<String>>) -> std::io::Result<()> {
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
                                var[0], 
                                var[1], 
                                var[2], 
                                var[3], 
                                var[4], 
                                var[5], 
                                var[6],
                                format!("{}{}",var[7],var[8]),
                                var[9],
                                var[10]);
        file.write(to_write.as_bytes()).expect("Error writing variant");
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use handlegraph::graph::HashGraph;
    use handlegraph::handlegraph::HandleGraph;
    use super::*;

    #[test]
    fn test_dfs_1() {
        let mut graph = HashGraph::new();
        let h1 = graph.append_handle("A");
        let h2 = graph.append_handle("CG");

        graph.create_edge(&h1, &h2);

        let mut g_dfs = HashGraph::new();
        dfs(&graph,&mut g_dfs,&h1.id());

        g_dfs.for_each_handle(|h| {
            display_node_edges(&g_dfs, &h);
            true
        });
        
    }

    #[test]
    fn test_dfs_2() {

        // Test a pattern that provides a different dfs_tree than Flavia's

        /*
        edges
        1  ----> 3 
          \-> 2 /
        */

        let mut graph = HashGraph::new();
        let h1 = graph.append_handle("A");
        let h2 = graph.append_handle("CG");
        let h3 = graph.append_handle("T");


        graph.create_edge(&h1, &h2);
        graph.create_edge(&h1, &h3);
        graph.create_edge(&h2, &h3);
        

        let mut g_dfs = HashGraph::new();
        dfs(&graph,&mut g_dfs,&h1.id());

        g_dfs.for_each_handle(|h| {
            display_node_edges(&g_dfs, &h);
            true
        });

        //println!("{:?}",graph);
        //println!("{:?}",g_dfs);
    }

    #[test]
    fn test_dfs_3() {

        // Test a pattern that provides the same dfs_tree as Flavia's

        /*
        edges
        
          /-> 3 -\ 
        1         4
          \-> 2 -/
        */

        let mut graph = HashGraph::new();
        let h1 = graph.append_handle("A");
        let h2 = graph.append_handle("CG");
        let h3 = graph.append_handle("T");
        let h4 = graph.append_handle("AC");

        graph.create_edge(&h1, &h2);
        graph.create_edge(&h1, &h3);
        graph.create_edge(&h2, &h4);
        graph.create_edge(&h3, &h4);

        let mut g_dfs = HashGraph::new();
        dfs(&graph,&mut g_dfs,&h1.id());

        //println!("{:?}",graph);

        g_dfs.for_each_handle(|h| {
            display_node_edges(&g_dfs, &h);
            true
        });

        //println!("{:?}",g_dfs);
    }
}