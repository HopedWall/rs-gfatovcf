use std::collections::HashMap;
use std::collections::BTreeMap; //like hashmap but sorted
use std::path::PathBuf;
use handlegraph::graph::HashGraph;
use handlegraph::graph::{PathId, PathStep};
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

//fn process_step(g : &HashGraph, s : &PathStep) -> String {
    //let h = g.get_handle_of_step(s);
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

//fn create_into_hashmap(g : &HashGraph, path_to_steps : &HashMap<String, Vec<String>>,  path : &PathId, step : &PathStep) -> bool {
    fn create_into_hashmap(g : &HashGraph, path_to_steps : &mut HashMap<String, Vec<String>>,  path : &PathId, step : &Handle) -> bool {    
    let path_name = g.get_path_name(path);
    println!("Print path name {}",path_name);
    if !path_to_steps.contains_key(path_name) {
        path_to_steps.insert(String::from(path_name), Vec::new());
    }
    println!("Something: {:?}",path_to_steps.get(path_name).unwrap());
    //.push(process_step(g, step));
    //path_to_steps[path_name].push(process_step(g,step));
    true
}


// TODO: fix last argument, currently removed
fn create_edge_and_so_on<F>(g : &HashGraph, g_dfs : &mut HashGraph, handle1 : &Handle, handle2 : &Handle, so_on_function : F) {
    let handle1_id = g.get_id(handle1);
    let handle2_id = g.get_id(handle2); 

    if !g_dfs.has_node(handle2_id) {
        //TODO: find a better solution than this
        //so_on_function(args[0],args[1]);
        
        if !g_dfs.has_node(handle1.id()) {
            g_dfs.create_handle(
                g.get_sequence(handle1), 
                handle1_id);
        }
    
        if !g_dfs.has_node(handle2.id()) {
            g_dfs.create_handle(
                g.get_sequence(handle2), 
                handle2_id);
        }
    
        g_dfs.create_edge(handle1, handle2)
    }
}
// //TODO: see create_edge_and_so_on
// fn dfs(g : &HashGraph, g_dfs : &mut HashGraph, node_id : &NodeId) {
//     let current_node = g.get_handle(*node_id, false);
//     let sequence_node = g.get_sequence(&current_node);

//     //TODO: fix this
//     g.follow_edges(
//         &current_node, 
//         Direction::Right,  //What should go here?
//         |neighbor| {
//             create_edge_and_so_on(
//                 &g, g_dfs, 
//                 &current_node, &neighbor, 
//                 dfs(g, g_dfs, &neighbor.id()));
//                 //&g.get_id(neighbor));
//                 true
//         });
// }

fn calculate_distance(visited_node_id_set : &mut HashSet<NodeId>, prev_node_id : &NodeId, neighbour_id : &NodeId, Q : &mut VecDeque<NodeId>, distances_map : &mut HashMap<NodeId,u64>) {
    if !visited_node_id_set.contains(&neighbour_id) {
        let previous_value = distances_map.get(prev_node_id).unwrap() + 1;
        distances_map.insert(*neighbour_id, previous_value);
        // is push_back correct?
        Q.push_back(*neighbour_id);
        visited_node_id_set.insert(*neighbour_id);
    }
}

// This function can be optimized via closures!
fn bfs_distances(g : &HashGraph, starting_node_id : &NodeId, distances_map : &HashMap<u64, u64>) -> (HashMap<NodeId,u64>, Vec<NodeId>) {
    let mut visited_node_id_set : HashSet<NodeId> = HashSet::new();
    let mut ordered_node_id_list = Vec::new();
    let mut distances_map : HashMap<NodeId, u64> = HashMap::new();
    let mut node_id_list = Vec::new();
    g.for_each_handle(|h| {
        node_id_list.push(g.get_id(h));
        true
    });
    for node_id in &node_id_list {
        distances_map.insert(*node_id, 0);
    }
    node_id_list.clear();

    let mut Q : VecDeque<NodeId> = VecDeque::new();
    Q.push_back(*starting_node_id);
    visited_node_id_set.insert(*starting_node_id);
    while !Q.is_empty() {
        let current_node_id = *Q.get(0).unwrap();
        let current_node = g.get_handle(current_node_id, false);
        ordered_node_id_list.push(current_node_id);

        g.follow_edges(
            &current_node, 
            Direction::Right, 
            |neighbor| {
                calculate_distance(&mut visited_node_id_set, &current_node_id, &g.get_id(neighbor), &mut Q, &mut distances_map);
                true
            });
    }

    (distances_map, ordered_node_id_list)
}

fn show_edge(g_dfs : &HashGraph, a : &Handle, b : &Handle) {
     print!("{} --> {}", g_dfs.get_id(a), g_dfs.get_id(b));
}

fn display_node_edges(g_dfs : &HashGraph, h : &Handle) {
    print!("node {}", g_dfs.get_id(h));
    g_dfs.follow_edges(
        h, Direction::Right,
        |n| {show_edge(&g_dfs,h,n); true}
    );
}

fn print_all_paths_util(g : &HashGraph, u : &NodeId, d : &NodeId, visited_node_id_set : &mut HashSet<NodeId>, path_list : &mut Vec<NodeId>, all_path_list : &mut Vec<NodeId>) {
    if !visited_node_id_set.contains(&u) {
        visited_node_id_set.insert(*u);
        path_list.push(*u);
    }

    //What does == mean exactly?
    if u == d {
        all_path_list.append(path_list);
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

fn print_all_paths(g : &HashGraph, start_node_id : &NodeId, end_node_id : &NodeId, all_path_list : &mut Vec<NodeId>) {
    let mut visited_node_id_set : HashSet<NodeId> = HashSet::new();
    
    let mut path_list : Vec<NodeId> = Vec::new();

    print_all_paths_util(g, start_node_id, end_node_id, &mut visited_node_id_set, &mut path_list, all_path_list);
}

fn main() {

    // let matches = App::new("rs-GFAtoVCF")
    //                       .version("1.0")
    //                       .author("Francesco Porto <francesco.porto97@gmail.com>")
    //                       .about("Converts GFA to VCF")
    //                       .arg(Arg::with_name("reference")
    //                            .short("r")
    //                            .long("reference")
    //                            .value_name("FILE")
    //                            .help("Sets the reference file to use")
    //                            .required(true)
    //                            .takes_value(true))
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
            println!("p: {}",p);
            graph.for_each_step_in_path(&p, 
                                        |s| {
                                            println!("p,s: {},{:?}",p,s);
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
        
    //     let mut node_id_to_path_and_pos_map : HashMap<NodeId, HashMap<String, usize>> = HashMap::new();
    //     for (path_name, steps_list) in &path_to_steps_map {
    //         let mut pos = 0;

    //         for nodeId_isRev in steps_list {
                
    //             // TODO: check what these are
    //             let node_id = nodeId_isRev.parse::<u64>().unwrap();
    //             let is_rev = nodeId_isRev;

    //             // Why does get handle require 2 parameters?
    //             let node_handle = graph.get_handle(NodeId::from(node_id), false);
    //             let seq = graph.get_sequence(&node_handle);

    //             let node_id_key = NodeId::from(node_id);

    //             if !node_id_to_path_and_pos_map.contains_key(&node_id_key) {
    //                 node_id_to_path_and_pos_map.insert(node_id_key, HashMap::new());
    //             }

    //             if !node_id_to_path_and_pos_map[&node_id_key].contains_key(path_name) {
    //                 //TODO: fix this
    //                 //(*node_id_to_path_and_pos_map.get(&node_id_key).unwrap()).insert(*path_name, pos);
    //                 //node_id_to_path_and_pos_map[&node_id_key][&path_name] = pos;
    //             }

    //             pos += seq.len();
    //         }
    //     }

    //     // node_id_to_path_and_pos_map must be sorted
    //     // TODO: find a better inplace solution
    //     let mut sorted = BTreeMap::new();
    //     for (id, value) in node_id_to_path_and_pos_map.iter() {
    //         sorted.insert(id, value);
    //     }
    //     for node_id in sorted.keys() {
    //         let path_and_pos_map = node_id_to_path_and_pos_map.get(node_id); 
    //         println!("Node_id : {}", node_id);

    //         //for (path, pos) in path_and_pos_map.into_iter() {
    //         //    println!("Path: {}  -- Pos: {}", path, pos);
    //         //}
    //     }

    //     let start_node = graph.get_handle(NodeId::from(1), false);

    //     let mut g_dfs = HashGraph::new();

    //     //Replace function, now in main directly
    //     //dfs(&graph, &mut g_dfs, &NodeId::from(1));
    //     let current_node = graph.get_handle(NodeId::from(1), false);
    //     let sequence_node = graph.get_sequence(&current_node);
    
    //     //TODO: fix this
    //     graph.follow_edges(
    //         &current_node, 
    //         Direction::Right,  //What should go here?
    //         |neighbor| {
    //             create_edge_and_so_on(
    //                 &graph, &mut g_dfs, 
    //                 &current_node, &neighbor, 
    //                 dfs(g, g_dfs, &neighbor.id()));
    //                 //&g.get_id(neighbor));
    //                 true
    //         });
        

    //     g_dfs.for_each_handle(|h| {
    //                                 display_node_edges(&g_dfs, &h);
    //                                 true
    //                               });
    //     let value = bfs_distances(&g_dfs, &NodeId::from(1), &HashMap::new()); 

    //     let distances_map = value.0;
    //     let ordered_node_id_list = value.1;

    //     for (node_id, distance) in distances_map.iter() {
    //         println!("{} - distance from root: {}", node_id, distance);
    //     }

    //     let mut dist_to_num_nodes : HashMap<u64,usize> = HashMap::new();

    //     for (_, distance) in distances_map.iter() {
    //         if !dist_to_num_nodes.contains_key(&distance) {
    //             dist_to_num_nodes.insert(*distance,0);
    //         }
    //         *dist_to_num_nodes.get_mut(distance).unwrap() += 1;
    //     }

    //     print!("Distance from root --> Num. nodes");
    //     for (k,v) in dist_to_num_nodes.iter() {
    //         println!("{} --> {}", k, v);
    //     }

    //     println!("Bubbles");
    //     let mut possible_bubbles_list : Vec<(NodeId,NodeId)> = Vec::new();
    //     let mut first_bubble = true;

    //     for node_id in ordered_node_id_list {
    //         let mut pair : (NodeId,NodeId) = (NodeId::from(0),NodeId::from(0));
    //         let key = distances_map[&node_id];
    //         if dist_to_num_nodes[&key] == 1 {
    //             if !first_bubble {
    //                 println!("{} END {:?} {}", node_id, node_id_to_path_and_pos_map[&node_id],g_dfs.get_sequence(&g_dfs.get_handle(node_id, false)));
                    
    //                 //possible_bubbles_list[possible_bubbles_list.len()][1] = node_id
    //                 pair.0 = node_id;
    //                 //possible_bubbles_list.push(node_id);
    //             }
    //             first_bubble = false;
    //             println!("{} START {:?} {}",node_id, node_id_to_path_and_pos_map[&node_id],g_dfs.get_sequence(&g_dfs.get_handle(node_id, false)));
    //             pair.1 = node_id;
    //             possible_bubbles_list.push(pair); //TODO: fix this
    //         } else {
    //             print!("{} Bolla {:?} {:?}",node_id, node_id_to_path_and_pos_map[&node_id], g_dfs.get_sequence(&g_dfs.get_handle(node_id,false)));
    //         }
    //     }

    //     println!("\n------------------");

    //     //TODO: fix this
    //     let mut path_to_sequence_map : HashMap<String,String> = HashMap::new();
    //     for (path_name, steps_list) in &path_to_steps_map {
    //         path_to_sequence_map.insert(path_name.to_string(), String::new());

    //          for node_id_rev in steps_list {
    //              //path_to_sequence_map[path_name].push_str(graph.get_sequence(&graph.get_handle(NodeId::from(node_id_rev.parse::<u64>().unwrap()), false)));
    //          }
    //    }

    //     let mut stuff_to_alts_dict : HashMap<String, HashSet<String>>  = HashMap::new();
    //     for current_ref in path_to_steps_map.keys() {
            
    //         let ref_path = Vec::new(); //TODO: fix this variable

    //         let length = 0;
    //         //let length = possible_bubbles_list.len();
    //         //TODO: needs to stop at -1
    //         //TODO: check this
    //         for (start,end) in &possible_bubbles_list {
                
    //             println!("ref_path: {:?}",ref_path);
    //             println!("Bubble [{},{}]",start, end);
                
    //             //TODO: Fix this
    //             //Find position of start in path
    //             //let start_node_index_in_ref_path = ref_path.iter().position(|&r| r == start).unwrap();
    //             let mut start_node_index_in_ref_path = 0;
    //             let mut all_path_list : Vec<Vec<NodeId>> = Vec::new();
    //             //print_all_paths(g, start, end, all_path_list)

    //             for path in all_path_list {
    //                 println!("\tPath: {:?}", path);
    //                 let mut pos_ref = node_id_to_path_and_pos_map[&start][current_ref]+1;
    //                 let mut pos_path = pos_ref;

    //                 println!("Start paths position: {}",pos_ref);
    //                 //TODO: fix this
    //                 //let max_index = cmp::min(path.len(), ref_path.len());
    //                 let mut max_index = ref_path.len();
    //                 let mut current_index_step_path = 0;
    //                 let mut current_index_step_ref = 0;

    //                 for i in 0..max_index {

    //                     let mut current_node_id_path = path[current_index_step_path];
    //                     let mut current_node_id_ref = ref_path[current_index_step_ref + start_node_index_in_ref_path];
                    
    //                     println!("{} {} ---> {} {}", pos_ref, pos_path, current_node_id_ref, current_node_id_path);
                        
    //                     if current_node_id_ref == current_node_id_path {
    //                         println!("REFERENCE");
    //                         //let node_seq = graph.get_sequence(&graph.get_handle(*current_node_id_ref, false));
    //                         //pos_ref += node_seq.len();
    //                         pos_path = pos_ref;
    //                         current_index_step_ref += 1;
    //                         current_index_step_path += 1;
    //                     } else {
    //                         let succ_node_id_path = path[current_index_step_path + 1];
    //                         let succ_node_id_ref = ref_path[current_index_step_ref + start_node_index_in_ref_path + 1];
    //                         if succ_node_id_ref == current_node_id_path {
    //                             println!("DEL");
    //                             let node_seq_ref = graph.get_sequence(&graph.get_handle(current_node_id_ref, false));
    //                             let prec_node_id_ref = ref_path[current_index_step_ref + start_node_index_in_ref_path - 1];
    //                             let prec_nod_seq_ref = graph.get_sequence(&graph.get_handle(prec_node_id_ref, false));
    //                             let key = [current_ref.to_string(), (pos_path - 1).to_string(), String::from(prec_nod_seq_ref), node_seq_ref.to_string()].join("_");
    //                             if !stuff_to_alts_dict.contains_key(&key) {
    //                                 stuff_to_alts_dict.insert(key, HashSet::new());
    //                             }
    //                             //stuff_to_alts_dict[key].add(prec_nod_seq_ref. + "_del");
    //                             pos_ref += node_seq_ref.len();
    //                             current_index_step_ref += 1;
    //                             current_node_id_ref = ref_path[current_index_step_ref + start_node_index_in_ref_path -1];
    //                             println!("\t {}", current_node_id_ref);
    //                             continue;
    //                         } else if succ_node_id_path == current_node_id_ref {
    //                             println!("INS");
    //                             let node_seq_path = graph.get_sequence(&graph.get_handle(current_node_id_path,false));
                                
    //                             let prec_node_id_ref = ref_path[current_index_step_ref + start_node_index_in_ref_path-1];
    //                             let prec_nod_seq_ref = graph.get_sequence(&graph.get_handle(prec_node_id_ref,false));
    //                             let key = [current_ref.to_string(), (pos_ref-1).to_string(), String::from(prec_nod_seq_ref)].join("-");
                                
    //                             if !stuff_to_alts_dict.contains_key(&key) {
    //                                 stuff_to_alts_dict.insert(key, HashSet::new());
    //                             }
    //                             //stuff_to_alts_dict[key].add(prec_nod_seq_ref[-1] + node_seq_path + '_ins')
    //                             pos_path += node_seq_path.len();
    //                             current_index_step_path += 1;
    //                             current_node_id_path = path[current_index_step_path];
    //                             println!("\t{}", current_node_id_path);
    //                             continue;
    //                         } else {
    //                             let node_seq_ref = graph.get_sequence(&graph.get_handle(current_node_id_ref, false));
    //                             let node_seq_path = graph.get_sequence(&graph.get_handle(current_node_id_path, false));

    //                             if node_seq_ref == node_seq_path {
    //                                 println!("REFERENCE");
    //                             } else {
    //                                 println!("SNV");
    //                             }

    //                             let key = [current_ref.to_string(), pos_path.to_string(), node_seq_ref.to_string()].join("-");
    //                             if !stuff_to_alts_dict.contains_key(&key){
    //                                 stuff_to_alts_dict.insert(key, HashSet::new());
    //                             }
                                
    //                             //stuff_to_alts_dict[&key].insert(String::from(node_seq_path).push_str("snv"));

    //                             pos_ref += node_seq_ref.len();
    //                             pos_path += node_seq_path.len();
    //                             current_index_step_ref += 1;
    //                             current_index_step_path += 1;
    //                         }
    //                     }
    //                 }
    //                 println!("---");
    //             }   
    //         }
    //         println!("==========================================");
    //    }


        //Find variations
        //let variations = find_variants(&graph, &reference);

        //Write variations to file
        //write_to_file(&PathBuf::from(out_path), &variations).unwrap();
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
        &["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"].join("\t"),
    ].join("\n");
    file.write(header.as_bytes()).expect("Error writing header");

    for var in variations {
        let to_write = format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                                var.chromosome, 
                                var.position, 
                                var.id, 
                                var.reference, 
                                var.alternate, 
                                var.quality, 
                                var.filter, 
                                var.info);
        file.write(to_write.as_bytes()).expect("Error writing variant");
    }
    Ok(())
}
