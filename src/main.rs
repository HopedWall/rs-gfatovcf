//! # GFAtoVCF
//! `GFAtoVCF` is a tool that finds variants in a Variation Graph.

use handlegraph::handle::{Direction, Edge, Handle, NodeId};
use handlegraph::handlegraph::HandleGraph;
use handlegraph::handlegraph::{handle_edges_iter, handles_iter};
use handlegraph::pathgraph::{paths_iter};
use handlegraph::hashgraph::HashGraph;
use handlegraph::hashgraph::PathId;
use handlegraph::mutablehandlegraph::MutableHandleGraph;
use handlegraph::pathgraph::steps_iter;
use handlegraph::pathgraph::PathHandleGraph;
use std::collections::BTreeMap; //like hashmap but sorted
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;
extern crate clap;
use clap::{App, Arg};
extern crate chrono;
use chrono::Utc;
use std::cmp;
use std::cmp::Ordering;
use std::collections::HashSet;
use std::collections::VecDeque;

use std::thread;
use std::sync::{Arc, Mutex};
use std::thread::JoinHandle;
use rayon::prelude::*;

/// A struct that holds Variants, as defined in the VCF format
#[derive(PartialEq)]
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
    sample_name: String,
}

/// Returns a step as a String with NodeId and Orientation
fn process_step(h: &Handle) -> String {
    let orient = if h.is_reverse() {
        "+".to_string()
    } else {
        "-".to_string()
    };

    format!("{}{}", h.id().to_string(), orient)
}

/// Returns all steps of a given path
fn get_steps_of_path(graph: &HashGraph, path_id : &PathId) -> Vec<String> {
    let mut steps_of_path : Vec<String> = vec![];
    graph.paths.get(path_id).unwrap().nodes.iter()
        .for_each(|handle| {
            steps_of_path.push(process_step(handle));
        });

    steps_of_path
}

/// Converts paths into sequences of nodes
fn paths_to_steps(graph : &HashGraph) -> HashMap<String, Vec<String>> {
    
    let mut path_to_steps_map: HashMap<String, Vec<String>> = HashMap::new();

    path_to_steps_map.par_extend(
        graph.paths.par_iter()
                            .map(|(path_id, _)| {
                                let path_name = graph.get_path(path_id).unwrap().name.clone();
                                let steps_list = get_steps_of_path(graph, &path_id);
                                (path_name,steps_list)
                            })
    );

    path_to_steps_map
}

// Wrapper function for bfs_new
fn bfs(g: &HashGraph, node_id: &NodeId) -> HashGraph {
    let mut g_bfs = HashGraph::new();
    bfs_support(g, &mut g_bfs, node_id);
    g_bfs
}

/// Computes the bfs of a given variation graph
fn bfs_support(g: &HashGraph, g_bfs: &mut HashGraph, node_id: &NodeId) {
    // let current_handle = Handle::pack(*node_id, false);
    // let mut added_handles: Vec<Handle> = vec![];

    // if !g_bfs.has_node(*node_id) {
    //     g_bfs.create_handle(g.sequence(current_handle), *node_id);
    // }

    // // Get neighbors for each node
    // let neighbors : Vec<_> = handle_edges_iter(g, current_handle, Direction::Right).collect();

    // neighbors.par_iter()
    //          .for_each(|neighbor| {
    //             bfs_parallel_step(g, g_bfs, &neighbor, &current_handle, &mut added_handles);
    //          });

    // //Repeat in parallel for all added handles
    // added_handles.par_iter()
    //             .for_each(|handle|{
    //                 bfs_support(g, g_bfs, &handle.id());
    //             });
        
}

fn bfs_parallel_step(g : &HashGraph, 
                    g_bfs : &mut HashGraph, 
                    neighbor : &Handle, 
                    current_handle : &Handle,
                    added_handles : &mut Vec<Handle>) {
    
    if !g_bfs.has_node(neighbor.id()) {
        let h = g_bfs.create_handle(g.sequence(*neighbor), neighbor.id());
        
        let edge = Edge::edge_handle(*current_handle, *neighbor);
        g_bfs.create_edge(&edge);
        
        //maybe better to return the edge to avoid deadlocks?
        added_handles.push(h);
    }
}

/// Prints an edge of a given HashGraph
fn show_edge(a: &Handle, b: &Handle) {
    println!("{} --> {}", a.id(), b.id());
}
/// Prints all nodes and edges of a given HashGraph
fn display_node_edges(g_dfs: &HashGraph, h: &Handle) {
    println!("node {}", h.id());

    //Obtain list of edges
    let edges : Vec<_> = handle_edges_iter(g_dfs, *h, Direction::Right).collect();
    
    edges.par_iter()    
        .for_each(|n| {
            show_edge(h, n)
        });   
    
}

/// Finds the distance of each node from a given root
fn bfs_distances(g: &HashGraph, starting_node_id: &NodeId) -> (BTreeMap<NodeId, u64>, Vec<NodeId>) {
    let mut visited_node_id_set = HashSet::new();
    let mut ordered_node_id_list = Vec::new();

    let mut distances_map = BTreeMap::new();
    let mut node_id_list = Vec::new();

    for handle in handles_iter(g) {
        let id = handle.id();
        node_id_list.push(id);
        distances_map.insert(id, 0);
    }

    let mut q: VecDeque<NodeId> = VecDeque::new();
    q.push_back(*starting_node_id);
    visited_node_id_set.insert(*starting_node_id);
    while let Some(current_node_id) = q.pop_front() {
        //let current_node_id = q.pop_front().unwrap();
        let current_node = Handle::pack(current_node_id, false);
        ordered_node_id_list.push(current_node_id);

        for h in handle_edges_iter(g, current_node, Direction::Right) {
            let n = h.id();
            if !visited_node_id_set.contains(&n) {
                let prev = *distances_map.get(&current_node_id).unwrap();
                distances_map.insert(n, prev + 1);
                q.push_back(n);
                visited_node_id_set.insert(n);
            }
        }
    }

    (distances_map, ordered_node_id_list)
}

/// Returns how many nodes are at the same distance
fn get_dist_to_num_nodes(distances_map: &BTreeMap<NodeId, u64>) -> BTreeMap<u64, usize> {
    let mut dist_to_num_nodes: BTreeMap<u64, usize> = BTreeMap::new();

    distances_map.values().for_each(|dist| {
        *dist_to_num_nodes.entry(*dist).or_default() += 1;
    });

    dist_to_num_nodes
}

/// Returns paths as sequences
fn get_path_to_sequence(
    graph: &HashGraph,
    path_to_steps_map: &HashMap<String, Vec<String>>,
) -> HashMap<String, String> {
    let mut path_to_sequence_map: HashMap<String, String> = HashMap::new();

    for (path_name, steps_list) in path_to_steps_map {
        path_to_sequence_map.insert(path_name.to_string(), String::new());

        for node_id_rev in steps_list {
            path_to_sequence_map
                .get_mut(path_name)
                .unwrap()
                .push_str(graph.sequence(Handle::pack(
                    NodeId::from(node_id_rev.parse::<u64>().unwrap()),
                    false,
                )));
        }
    }

    path_to_sequence_map
}

/// Detects the bubbles in a variation graph
fn detect_bubbles(
    distances_map: &BTreeMap<NodeId, u64>,
    ordered_node_id_list: &[NodeId],
    dist_to_num_nodes: &BTreeMap<u64, usize>,
) -> Vec<(NodeId, NodeId)> {
    let mut possible_bubbles_list: Vec<(NodeId, NodeId)> = Vec::new();
    let mut first_bubble = true;

    for node_id in ordered_node_id_list {
        let mut pair: (NodeId, NodeId) = (NodeId::from(0), NodeId::from(0));
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

    //Delete last bubble, won't be used anyways
    possible_bubbles_list.pop();

    possible_bubbles_list
}

/// Prints all paths between two nodes
fn print_all_paths_util(
    g: &HashGraph,
    u: &NodeId,
    d: &NodeId,
    visited_node_id_set: &mut HashSet<NodeId>,
    path_list: &mut Vec<NodeId>,
    all_path_list: &mut Vec<Vec<NodeId>>,
) {
    if !visited_node_id_set.contains(&u) {
        visited_node_id_set.insert(*u);
        path_list.push(*u);
    }

    if u == d {
        all_path_list.push(path_list.to_vec());
    } else {
        for i_node in handle_edges_iter(g, Handle::pack(*u, false), Direction::Right) {
            print_all_paths_util(
                g,
                &i_node.id(),
                d,
                visited_node_id_set,
                path_list,
                all_path_list,
            );
        }
    }

    path_list.pop();
    visited_node_id_set.remove(&u);
}

/// Prints all paths in a given HashGrap, starting from a specific node and ending in another node
fn print_all_paths(
    g: &HashGraph,
    start_node_id: &NodeId,
    end_node_id: &NodeId,
    all_path_list: &mut Vec<Vec<NodeId>>,
) {
    let mut visited_node_id_set: HashSet<NodeId> = HashSet::new();

    let mut path_list: Vec<NodeId> = Vec::new();

    // let start_handle = Handle::pack(*start_node_id, false);
    // let end_handle = Handle::pack(*end_node_id, false);

    // for neighbor in handle_edges_iter(g, start_handle, Direction::Right) {

    //     if !visited_node_id_set.contains(&neighbor.id()) {
    //         visited_node_id_set.insert(neighbor.id());
    //         path_list.push(neighbor.id());
    //     }

    //     if neighbor == end_handle {
    //         all_path_list.push(path_list.to_vec());
    //     } else {
    //         for n_of_n in handle_edges_iter(g, neighbor, Direction::Right) {

    //         }
    //     }
    // }

    print_all_paths_util(
        g,
        start_node_id,
        end_node_id,
        &mut visited_node_id_set,
        &mut path_list,
        all_path_list,
    );
}

/// Detects variants in a variation graph
fn detect_all_variants(
    path_to_steps_map: &HashMap<String, Vec<String>>,
    possible_bubbles_list: &[(NodeId, NodeId)],
    graph: &HashGraph,
    node_id_to_path_and_pos_map: &BTreeMap<NodeId, HashMap<String, usize>>,
    verbose: bool,
) -> Vec<Variant> {
    let mut stuff_to_alts_map: HashMap<String, HashSet<String>> = HashMap::new();

    // Taking each known path as reference, explore all bubbles in order to find variants;
    // these will be stored in stuff_to_alts_map
    for current_ref in path_to_steps_map.keys() {
        //Obtain all steps for each path
        //let mut ref_path = vec![];
        let ref_path: Vec<u64> = path_to_steps_map[current_ref]
            .iter()
            .map(|x| x.parse::<u64>().unwrap())
            .collect();

        if verbose {
            println!("path_to_steps_map: {:?}", path_to_steps_map);
        }

        // for x in &path_to_steps_map[current_ref] {
        //     ref_path.push(x.parse::<u64>().unwrap());
        // }

        detect_variants_per_reference(
            &current_ref,
            &ref_path,
            possible_bubbles_list,
            graph,
            node_id_to_path_and_pos_map,
            &mut stuff_to_alts_map,
            verbose,
        );
    }

    // Convert stuff_to_alts_map to a more readable format
    let mut vcf_list: Vec<Variant> = Vec::new();
    for (chrom_pos_ref, alt_type_set) in &stuff_to_alts_map {
        let vec: Vec<&str> = chrom_pos_ref.split('_').collect();
        let chrom = vec[0];
        let pos = vec[1];
        let refr = vec[2];

        let (alt_list, type_set): (Vec<_>, Vec<_>) = alt_type_set
            .iter()
            .map(|x| {
                let split: Vec<_> = x.split('_').collect();
                (split[0].to_string(), split[1])
            })
            .unzip();

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
    vcf_list.sort_by(|a, b| match a.chromosome.cmp(&b.chromosome) {
        Ordering::Equal => a
            .position
            .parse::<i32>()
            .unwrap()
            .cmp(&b.position.parse::<i32>().unwrap()),
        other => other,
    });

    vcf_list
}
/// Detect variants for a specific reference
fn detect_variants_per_reference(
    current_ref: &str,
    ref_path: &[u64],
    possible_bubbles_list: &[(NodeId, NodeId)],
    graph: &HashGraph,
    node_id_to_path_and_pos_map: &BTreeMap<NodeId, HashMap<String, usize>>,
    stuff_to_alts_map: &mut HashMap<String, HashSet<String>>,
    verbose: bool,
) {
    // Create closure that will be used later
    let get_last = |prec_node_seq_ref: &str, node_seq_ref| {
        let mut last = prec_node_seq_ref[prec_node_seq_ref.len() - 1..].to_string();
        last.push_str(node_seq_ref);
        last
    };

    // Check all bubbles
    for (start, end) in possible_bubbles_list {
        if verbose {
            println!("ref_path: {:?}", ref_path);
            println!("Bubble [{},{}]", start, end);
        }

        let start_node_index_in_ref_path: usize;
        match ref_path.iter().position(|&r| NodeId::from(r) == *start) {
            None => continue, //ignore, start not found in ref path
            Some(r) => start_node_index_in_ref_path = r,
        };

        let mut all_path_list: Vec<Vec<NodeId>> = Vec::new();
        print_all_paths(&graph, start, end, &mut all_path_list);

        //println!("All paths list: {:?}",all_path_list);
        for path in &all_path_list {
            if verbose {
                println!("\tPath: {:?}", path);
            }

            let mut pos_ref = node_id_to_path_and_pos_map[start][current_ref] + 1;
            let mut pos_path = pos_ref;

            let max_index = cmp::min(path.len(), ref_path.len());

            let mut current_index_step_path = 0;
            let mut current_index_step_ref = 0;

            for _i in 0..max_index {
                let mut current_node_id_ref =
                    NodeId::from(ref_path[current_index_step_ref + start_node_index_in_ref_path]);
                let mut current_node_id_path = path[current_index_step_path];

                if verbose {
                    println!(
                        "{} {} ---> {} {}",
                        pos_ref, pos_path, current_node_id_ref, current_node_id_path
                    );
                }

                if current_node_id_ref == current_node_id_path {
                    if verbose {
                        println!("REFERENCE");
                    }

                    let node_seq = graph.sequence(Handle::pack(current_node_id_ref, false));
                    pos_ref += node_seq.len();
                    pos_path = pos_ref;

                    current_index_step_ref += 1;
                    current_index_step_path += 1;
                } else {
                    // Shouldn't be happening anymore
                    if current_index_step_path + 1 >= path.len() {
                        break;
                    }
                    if current_index_step_ref + start_node_index_in_ref_path + 1 >= ref_path.len() {
                        break;
                    }

                    let succ_node_id_path = path[current_index_step_path + 1];
                    let succ_node_id_ref = NodeId::from(
                        ref_path[current_index_step_ref + start_node_index_in_ref_path + 1],
                    );
                    if succ_node_id_ref == current_node_id_path {
                        if verbose {
                            println!("DEL");
                        }

                        let node_seq_ref = graph.sequence(Handle::pack(current_node_id_ref, false));

                        let prec_node_id_ref = NodeId::from(
                            ref_path[current_index_step_ref + start_node_index_in_ref_path - 1],
                        );
                        let prec_nod_seq_ref =
                            graph.sequence(Handle::pack(prec_node_id_ref, false));

                        let last = get_last(prec_nod_seq_ref, node_seq_ref);

                        let key =
                            [current_ref.to_string(), (pos_path - 1).to_string(), last].join("_");
                        stuff_to_alts_map.entry(key).or_insert(HashSet::new());
                        //TODO: find a better way to do this
                        let last = get_last(prec_nod_seq_ref, node_seq_ref);

                        let key =
                            [current_ref.to_string(), (pos_path - 1).to_string(), last].join("_");

                        let last = prec_nod_seq_ref[prec_nod_seq_ref.len() - 1..].to_string();
                        let mut string_to_insert = last;
                        string_to_insert.push_str("_del");
                        stuff_to_alts_map
                            .get_mut(&key)
                            .unwrap()
                            .insert(string_to_insert);

                        pos_ref += node_seq_ref.len();

                        current_index_step_ref += 1;
                        current_node_id_ref = NodeId::from(
                            ref_path[current_index_step_ref + start_node_index_in_ref_path - 1],
                        );
                        if verbose {
                            println!("\t {}", current_node_id_ref);
                        }

                        continue;
                    } else if succ_node_id_path == current_node_id_ref {
                        if verbose {
                            println!("INS");
                        }

                        let node_seq_path =
                            graph.sequence(Handle::pack(current_node_id_path, false));

                        let prec_node_id_ref = NodeId::from(
                            ref_path[current_index_step_ref + start_node_index_in_ref_path - 1],
                        );
                        let prec_nod_seq_ref =
                            graph.sequence(Handle::pack(prec_node_id_ref, false));

                        let last = prec_nod_seq_ref[prec_nod_seq_ref.len() - 1..].to_string();
                        //let key = [current_ref.to_string(), (pos_ref-1).to_string(), String::from(prec_nod_seq_ref)].join("_");
                        let key =
                            [current_ref.to_string(), (pos_ref - 1).to_string(), last].join("_");

                        stuff_to_alts_map.entry(key).or_insert(HashSet::new());

                        //Re-create key since it goes out of scope
                        let last = prec_nod_seq_ref[prec_nod_seq_ref.len() - 1..].to_string();
                        let key =
                            [current_ref.to_string(), (pos_ref - 1).to_string(), last].join("_");

                        let last = prec_nod_seq_ref[prec_nod_seq_ref.len() - 1..].to_string();
                        let mut string_to_insert = last;
                        string_to_insert.push_str(node_seq_path);
                        string_to_insert.push_str("_ins");
                        stuff_to_alts_map
                            .get_mut(&key)
                            .unwrap()
                            .insert(string_to_insert);

                        pos_path += node_seq_path.len();

                        current_index_step_path += 1;
                        current_node_id_path = path[current_index_step_path];
                        if verbose {
                            println!("\t{}", current_node_id_path);
                        }

                        continue;
                    } else {
                        let node_seq_ref = graph.sequence(Handle::pack(current_node_id_ref, false));
                        let node_seq_path =
                            graph.sequence(Handle::pack(current_node_id_path, false));

                        if node_seq_ref == node_seq_path {
                            if verbose {
                                println!("REFERENCE");
                            }
                        } else {
                            if verbose {
                                println!("SNV");
                            }
                        }

                        let key = [
                            current_ref.to_string(),
                            pos_path.to_string(),
                            node_seq_ref.to_string(),
                        ]
                        .join("_");

                        stuff_to_alts_map.entry(key).or_insert(HashSet::new());

                        //TODO: find a better way to do this
                        let key = [
                            current_ref.to_string(),
                            pos_path.to_string(),
                            node_seq_ref.to_string(),
                        ]
                        .join("_");
                        let mut string_to_insert =
                            node_seq_path.chars().last().unwrap().to_string();
                        string_to_insert.push_str("_snv");
                        stuff_to_alts_map
                            .get_mut(&key)
                            .unwrap()
                            .insert(string_to_insert);

                        pos_ref += node_seq_ref.len();
                        pos_path += node_seq_path.len();
                        current_index_step_ref += 1;
                        current_index_step_path += 1;
                    }
                }
            }
            if verbose {
                println!("---");
            }
        }
    }
    if verbose {
        println!("==========================================");
    }
}

/// Gets the position of a node in all the paths it is present in
fn get_node_positions_in_paths(
    graph: &HashGraph,
    path_to_steps_map: &mut HashMap<String, Vec<String>>,
) -> BTreeMap<NodeId, HashMap<String, usize>> {
    let mut node_id_to_path_and_pos_map: BTreeMap<NodeId, HashMap<String, usize>> = BTreeMap::new();
    //let mut b: BTreeMap<NodeId, HashMap<String, usize>> = BTreeMap::new();

    node_id_to_path_and_pos_map.par_extend(
       //b.into_par_iter()
        
        path_to_steps_map.par_iter()
              .for_each(|(path_name, steps_list)| {
              
                let mut pos = 0;

                steps_list.iter_mut()
                    .for_each(|node_id_is_rev|{
                        
                        // Get orientation
                        let _is_rev = node_id_is_rev.pop().unwrap();
                        
                        // Get the id of the node string -> NodeId
                        let node_id: NodeId = NodeId::from(node_id_is_rev.parse::<u64>().unwrap());

                        let node_handle = Handle::pack(node_id, false);
                        let seq = graph.sequence(node_handle);

                        let mut n = node_id_to_path_and_pos_map;
                        n.entry(node_id).or_insert(HashMap::new());

                        if !n[&node_id].contains_key(path_name) {
                                n
                                .get_mut(&node_id)
                                .unwrap()
                                .insert(String::from(path_name), pos);
                        }

                        pos += seq.len();
                    
                    });
                    
        })
    );

    node_id_to_path_and_pos_map
}

/// The function that runs the script
fn main() {
    let matches = App::new("rs-GFAtoVCF")
        .version("1.0")
        .author("Francesco Porto <francesco.porto97@gmail.com>")
        .about("Converts GFA to VCF")
        .arg(
            Arg::with_name("input")
                .short("i")
                .long("input")
                .value_name("FILE")
                .help("Sets the input file to use")
                .required(true)
                .takes_value(true),
        )
        .arg(
            Arg::with_name("output")
                .short("o")
                .long("output")
                .value_name("FILE")
                .help("Sets the output file to use")
                .required(true)
                .takes_value(true),
        )
        .arg(
            Arg::with_name("verbose")
                .short("v")
                .long("verbose")
                .help("Sets whether to display debug messages or not"),
        )
        .get_matches();

    let in_path_file = matches
        .value_of("input")
        .expect("Could not parse argument --input");
    let out_path_file = matches
        .value_of("output")
        .expect("Could not parse argument --output");
    let verbose = matches.is_present("verbose");

    //let in_path_file = "./input/samplePath3.gfa";
    //let out_path_file = "./input/samplePath3.vcf";

    if let Some(gfa) = gfa::parser::parse_gfa(&PathBuf::from(in_path_file)) {
        let graph = HashGraph::from_gfa(&gfa);

        //Obtains, for each path, a list of all its steps (its nodes)
        //let path_to_steps_map: HashMap<String, Vec<String>> = paths_to_steps(&graph);
        let mut path_to_steps_map = paths_to_steps(&graph);
        println!("Path to steps map: {:#?}",path_to_steps_map);

        // Obtains, for each node, its position in each path where the node is in
        let node_id_to_path_and_pos_map: BTreeMap<NodeId, HashMap<String, usize>> =
              get_node_positions_in_paths(&graph, &mut path_to_steps_map);

        // if verbose {
        //     for node_id in node_id_to_path_and_pos_map.keys() {
        //         let path_and_pos_map = node_id_to_path_and_pos_map.get(node_id);
        //         println!("Node_id : {}", node_id);

        //         for (path, pos) in path_and_pos_map.unwrap() {
        //             println!("Path: {}  -- Pos: {}", path, pos);
        //         }
        //     }
        // }

        //Obtains the tree representing the bfs
        //let graph_arc = Arc::new(graph);
        //let g_bfs: HashGraph = bfs(&graph, &NodeId::from(1));
        

        // if verbose {
        //     handles_iter(graph_arc.as_ref())
        //         .for_each(|h| {
        //             let cloned = Arc::clone(&graph_arc);
        //             display_node_edges(cloned, Arc::new(h));
        //         });    
        // }

        // // Obtains, for each level of the tree, how many nodes are there
        // let (distances_map, ordered_node_id_list) = bfs_distances(&g_bfs, &NodeId::from(1));

        // if verbose {
        //     println!("\nNode --> Distance from root");
        //     for (node_id, distance) in distances_map.iter() {
        //         println!("{} - distance from root: {}", node_id, distance);
        //     }
        // }

        // //Obtains a map where, for each distance from root, the number of nodes
        // //at that distance are present
        // let dist_to_num_nodes: BTreeMap<u64, usize> = get_dist_to_num_nodes(&distances_map);

        // if verbose {
        //     println!("\nDistance from root --> Num. nodes");
        //     for (k, v) in dist_to_num_nodes.iter() {
        //         println!("{} --> {}", k, v);
        //     }
        // }

        // let possible_bubbles_list: Vec<(NodeId, NodeId)> =
        //     detect_bubbles(&distances_map, &ordered_node_id_list, &dist_to_num_nodes);

        // if verbose {
        //     println!("\nBubbles");
        //     println!("Detected bubbles {:#?}", possible_bubbles_list);
        //     println!("\n------------------");

        //     //Obtains, for each path, a string representing all the bases of the path (not actually used)
        //     let path_to_sequence_map: HashMap<String, String> =
        //         get_path_to_sequence(&graph, &path_to_steps_map);
        //     println!("Path to sequence: {:?}", path_to_sequence_map);
        // }

        // let vcf_list = detect_all_variants(
        //     &path_to_steps_map,
        //     &possible_bubbles_list,
        //     &graph,
        //     &node_id_to_path_and_pos_map,
        //     verbose,
        // );

        // //Write variants to file
        // write_to_file(&PathBuf::from(out_path_file), &vcf_list).unwrap();
    } else {
        panic!("Couldn't parse gfa file!");
    }
}

/// Write variants to file
fn write_to_file(path: &PathBuf, variants: &[Variant]) -> std::io::Result<()> {
    let mut file = File::create(path).unwrap_or_else(|_| panic!("Error creating file {:?}", path));

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
    file.write_all(header.as_bytes())
        .expect("Error writing header");

    file.write_all(b"\n").expect("Error writing to file");
    for var in variants {
        let to_write = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            var.chromosome,
            var.position,
            var.id,
            var.reference,
            var.alternate,
            var.quality,
            var.filter,
            var.info,
            var.format,
            var.sample_name
        );
        file.write_all(to_write.as_bytes())
            .expect("Error writing variant");
    }
    Ok(())
}

// #[cfg(test)]
// mod tests {
//     use super::*;
//     use handlegraph::handlegraph::HandleGraph;
//     use handlegraph::hashgraph::HashGraph;

//     //Used in other tests
//     fn read_test_gfa() -> HashGraph {
//         use gfa::parser::parse_gfa;

//         HashGraph::from_gfa(&parse_gfa(&PathBuf::from("./input/samplePath3.gfa")).unwrap())
//     }

//     #[test]
//     fn test_path_to_steps() {
//         let graph = read_test_gfa();

//         let path_to_steps_map: HashMap<String, Vec<String>> = paths_to_steps(&graph);

//         // Check if all paths have been found
//         assert_eq!(path_to_steps_map.keys().len(), 3);
//         assert!(path_to_steps_map.contains_key("x"));
//         assert!(path_to_steps_map.contains_key("y"));
//         assert!(path_to_steps_map.contains_key("z"));

//         //Check for each path that all its node have been found

//         //P	x	1+,3+,5+,6+,8+,9+,11+,12+,13+,15+,16+,18+,19+
//         let path_x: Vec<String> = vec![
//             "1+".to_string(),
//             "3+".to_string(),
//             "5+".to_string(),
//             "6+".to_string(),
//             "8+".to_string(),
//             "9+".to_string(),
//             "11+".to_string(),
//             "12+".to_string(),
//             "13+".to_string(),
//             "15+".to_string(),
//             "16+".to_string(),
//             "18+".to_string(),
//             "19+".to_string(),
//         ];
//         assert_eq!(*path_to_steps_map.get("x").unwrap(), path_x);

//         //P	y	1+,2+,5+,6+,8+,9+,10+,11+,13+,14+,16+,17+,19+
//         let path_y: Vec<String> = vec![
//             "1+".to_string(),
//             "2+".to_string(),
//             "5+".to_string(),
//             "6+".to_string(),
//             "8+".to_string(),
//             "9+".to_string(),
//             "10+".to_string(),
//             "11+".to_string(),
//             "13+".to_string(),
//             "14+".to_string(),
//             "16+".to_string(),
//             "17+".to_string(),
//             "19+".to_string(),
//         ];
//         assert_eq!(*path_to_steps_map.get("y").unwrap(), path_y);

//         //P	z	1+,2+,4+,6+,7+,9+,10+,11+,12+,13+,14+,16+,18+,19+
//         let path_z: Vec<String> = vec![
//             "1+".to_string(),
//             "2+".to_string(),
//             "4+".to_string(),
//             "6+".to_string(),
//             "7+".to_string(),
//             "9+".to_string(),
//             "10+".to_string(),
//             "11+".to_string(),
//             "12+".to_string(),
//             "13+".to_string(),
//             "14+".to_string(),
//             "16+".to_string(),
//             "18+".to_string(),
//             "19+".to_string(),
//         ];
//         assert_eq!(*path_to_steps_map.get("z").unwrap(), path_z);
//     }

//     #[test]
//     fn test_bfs() {
//         let graph = read_test_gfa();
//         let g_bfs = bfs(&graph, &NodeId::from(1));

//         // All nodes must be present in bfs
//         assert_eq!(graph.node_count(), g_bfs.node_count());

//         // There should be less (or the same number of) edges in g_bfs
//         // since nodes can only get added once
//         assert!(graph.edge_count() >= g_bfs.edge_count());
//     }

//     #[test]
//     fn test_bfs_distances() {
//         let graph = read_test_gfa();
//         let g_bfs = bfs(&graph, &NodeId::from(1));

//         let (distances_map, ordered_node_id_list) = bfs_distances(&g_bfs, &NodeId::from(1));

//         // Should pass easily
//         assert_eq!(
//             distances_map[&NodeId::from(2)],
//             distances_map[&NodeId::from(3)]
//         );
//         assert_eq!(
//             distances_map[&NodeId::from(4)],
//             distances_map[&NodeId::from(5)]
//         );
//         assert_eq!(
//             distances_map[&NodeId::from(7)],
//             distances_map[&NodeId::from(8)]
//         );

//         // Caused problems in dfs
//         assert_eq!(
//             distances_map[&NodeId::from(10)],
//             distances_map[&NodeId::from(11)]
//         );
//         assert_eq!(
//             distances_map[&NodeId::from(12)],
//             distances_map[&NodeId::from(13)]
//         );

//         //Should pass easily
//         assert_eq!(
//             distances_map[&NodeId::from(14)],
//             distances_map[&NodeId::from(15)]
//         );
//         assert_eq!(
//             distances_map[&NodeId::from(17)],
//             distances_map[&NodeId::from(18)]
//         );

//         // Check ordered_node_list now
//         let mut sorted = ordered_node_id_list.clone();
//         sorted.sort();
//         assert!(ordered_node_id_list.eq(&sorted));
//     }

//     #[test]
//     fn test_dist_to_num_nodes() {
//         let graph = read_test_gfa();
//         let g_bfs = bfs(&graph, &NodeId::from(1));

//         let (distances_map, _) = bfs_distances(&g_bfs, &NodeId::from(1));
//         let mut dist_to_num_nodes: BTreeMap<u64, usize> = BTreeMap::new();
//         for (_, distance) in distances_map.iter() {
//             if !dist_to_num_nodes.contains_key(&distance) {
//                 dist_to_num_nodes.insert(*distance, 0);
//             }
//             *dist_to_num_nodes.get_mut(distance).unwrap() += 1;
//         }

//         //println!("{:#?}",dist_to_num_nodes);

//         //Root
//         assert_eq!(dist_to_num_nodes[&0], 1);

//         assert_eq!(dist_to_num_nodes[&1], 2);
//         assert_eq!(dist_to_num_nodes[&2], 2);
//         assert_eq!(dist_to_num_nodes[&3], 1);
//         assert_eq!(dist_to_num_nodes[&4], 2);
//         assert_eq!(dist_to_num_nodes[&5], 1);

//         // Critical, did not work with dfs
//         assert_eq!(dist_to_num_nodes[&6], 2);
//         assert_eq!(dist_to_num_nodes[&7], 2);
//         assert_eq!(dist_to_num_nodes[&8], 2);

//         assert_eq!(dist_to_num_nodes[&9], 1);
//         assert_eq!(dist_to_num_nodes[&10], 2);
//         assert_eq!(dist_to_num_nodes[&11], 1);
//     }

//     #[test]
//     fn test_bubble_detection() {
//         let graph = read_test_gfa();
//         let g_bfs = bfs(&graph, &NodeId::from(1));

//         let (distances_map, ordered_node_id_list) = bfs_distances(&g_bfs, &NodeId::from(1));
//         let mut dist_to_num_nodes: BTreeMap<u64, usize> = BTreeMap::new();
//         for (_, distance) in distances_map.iter() {
//             if !dist_to_num_nodes.contains_key(&distance) {
//                 dist_to_num_nodes.insert(*distance, 0);
//             }
//             *dist_to_num_nodes.get_mut(distance).unwrap() += 1;
//         }

//         let possible_bubbles_list: Vec<(NodeId, NodeId)> =
//             detect_bubbles(&distances_map, &ordered_node_id_list, &dist_to_num_nodes);

//         //println!("Possible bubbles list: {:#?}",possible_bubbles_list);

//         assert!(possible_bubbles_list.contains(&(NodeId::from(1), NodeId::from(6))));
//         assert!(possible_bubbles_list.contains(&(NodeId::from(6), NodeId::from(9))));
//         assert!(possible_bubbles_list.contains(&(NodeId::from(9), NodeId::from(16))));
//         assert!(possible_bubbles_list.contains(&(NodeId::from(16), NodeId::from(19))));
//     }

//     fn run_whole_script() -> Vec<Variant> {
//         let graph = read_test_gfa();

//         //Obtain preliminary data required for future steps
//         let mut path_to_steps_map: HashMap<String, Vec<String>> = paths_to_steps(&graph);
//         let node_id_to_path_and_pos_map: BTreeMap<NodeId, HashMap<String, usize>> =
//             get_node_positions_in_paths(&graph, &mut path_to_steps_map);

//         //Compute bfs and analyze the results
//         let g_bfs: HashGraph = bfs(&graph, &NodeId::from(1));
//         let (distances_map, ordered_node_id_list) = bfs_distances(&g_bfs, &NodeId::from(1));
//         let dist_to_num_nodes: BTreeMap<u64, usize> = get_dist_to_num_nodes(&distances_map);

//         //Find the bubbles
//         let possible_bubbles_list: Vec<(NodeId, NodeId)> =
//             detect_bubbles(&distances_map, &ordered_node_id_list, &dist_to_num_nodes);

//         //Find variants from bubbles
//         let vcf_list = detect_all_variants(
//             &path_to_steps_map,
//             &possible_bubbles_list,
//             &graph,
//             &node_id_to_path_and_pos_map,
//             false,
//         );

//         vcf_list
//     }

//     #[test]
//     fn test_variant_detection() {
//         let variants_found = run_whole_script();

//         //Check that all variants have been found
//         assert_eq!(variants_found.len(), 21);

//         //Check one variant per type

//         //x	9	.	G	A	.	.	TYPE=snv	GT	0|1
//         let snv = Variant {
//             chromosome: "x".to_string(),
//             position: "9".to_string(),
//             id: ".".to_string(),
//             reference: "G".to_string(),
//             alternate: "A".to_string(),
//             quality: ".".to_string(),
//             filter: ".".to_string(),
//             info: format!("TYPE={}", "snv"),
//             format: "GT".to_string(),
//             sample_name: "0|1".to_string(),
//         };
//         assert!(variants_found.contains(&snv));

//         //x	18	.	T	TAA	.	.	TYPE=ins	GT	0|1
//         let ins = Variant {
//             chromosome: "x".to_string(),
//             position: "18".to_string(),
//             id: ".".to_string(),
//             reference: "T".to_string(),
//             alternate: "TAA".to_string(),
//             quality: ".".to_string(),
//             filter: ".".to_string(),
//             info: format!("TYPE={}", "ins"),
//             format: "GT".to_string(),
//             sample_name: "0|1".to_string(),
//         };
//         assert!(variants_found.contains(&ins));

//         //y	18	.	TAA	T	.	.	TYPE=del	GT	0|1
//         let del = Variant {
//             chromosome: "y".to_string(),
//             position: "18".to_string(),
//             id: ".".to_string(),
//             reference: "TAA".to_string(),
//             alternate: "T".to_string(),
//             quality: ".".to_string(),
//             filter: ".".to_string(),
//             info: format!("TYPE={}", "del"),
//             format: "GT".to_string(),
//             sample_name: "0|1".to_string(),
//         };
//         assert!(variants_found.contains(&del));
//     }
// }
