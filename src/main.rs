//! # GFAtoVCF
//! `GFAtoVCF` is a tool that finds variants in a Variation Graph.

use handlegraph::handle::{Direction, Edge, Handle, NodeId};
use handlegraph::handlegraph::HandleGraph;
use handlegraph::handlegraph::{handle_edges_iter, handles_iter};
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
//use std::io;
use json::*;
use std::io::BufWriter;

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
/// Returns all paths as a hashmap, having the path_name as key and a list of steps as values
fn create_into_hashmap(
    g: &HashGraph,
    path_to_steps: &mut HashMap<String, Vec<String>>,
    path: &PathId,
    step: &Handle,
) -> bool {
    let path_name = g.get_path(path).unwrap().name.clone();

    path_to_steps
        .entry(path_name)
        .or_default()
        .push(process_step(step));

    true
}
/// Converts paths into sequences of nodes
fn paths_to_steps(graph: &HashGraph) -> HashMap<String, Vec<String>> {
    let mut path_to_steps_map: HashMap<String, Vec<String>> = HashMap::new();

    for path_id in std::iter::from_fn(graph.paths_iter_impl()) {
        for step in steps_iter(graph, path_id) {
            let handle = graph.handle_of_step(&step).unwrap();
            create_into_hashmap(graph, &mut path_to_steps_map, path_id, &handle);
        }
    }

    path_to_steps_map
}

/// Wrapper function for bfs_new
fn bfs(g: &HashGraph, node_id: &NodeId) -> HashGraph {
    let mut g_bfs = HashGraph::new();

    //Create queue
    let mut q: VecDeque<NodeId> = VecDeque::new();

    //Insert first value
    q.push_back(*node_id);

    while !q.is_empty() {
        //println!("Queue is {:#?}",q);

        let curr_node = q.pop_front().unwrap();
        let current_handle = Handle::pack(curr_node, false);

        //Check if curr_node is already in g_bfs
        if !g_bfs.has_node(curr_node) {
            g_bfs.create_handle(g.sequence(current_handle), curr_node);
        }

        for neighbor in handle_edges_iter(g, current_handle, Direction::Right) {
            if !g_bfs.has_node(neighbor.id()) {
                //Create handle in g_bfs
                g_bfs.create_handle(g.sequence(neighbor), neighbor.id());

                //Add neighbor id to queue
                q.push_back(neighbor.id());

                //Create edge from curr_handle to new node in g_bfs
                let edge = Edge::edge_handle(current_handle, neighbor);
                g_bfs.create_edge(&edge);

                //Add new node to queue
                q.push_back(neighbor.id());
            }
        }
    }

    g_bfs
}

/// Prints an edge of a given HashGraph
fn show_edge(a: &Handle, b: &Handle) {
    println!("{} --> {}", a.id(), b.id());
}
/// Prints all nodes and edges of a given HashGraph
fn display_node_edges(g_dfs: &HashGraph, h: &Handle) {
    println!("node {}", h.id());

    for n in handle_edges_iter(g_dfs, *h, Direction::Right) {
        show_edge(h, &n);
    }
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

/// Prints all paths in a given HashGrap, starting from a specific node and ending in another node
fn find_all_paths_between(
    g: &HashGraph,
    start_node_id: &NodeId,
    end_node_id: &NodeId,
) -> Vec<Vec<NodeId>> {
    let mut all_paths_list: Vec<Vec<NodeId>> = Vec::new();

    let mut visited_node_id_set: HashSet<NodeId> = HashSet::new();
    //Create queue
    let mut q: VecDeque<NodeId> = VecDeque::new();
    //Insert first value
    q.push_back(*start_node_id);
    all_paths_list.push(vec![*start_node_id]);

    while !q.is_empty() {
        //println!("All paths is {:#?}",all_paths_list);
        //println!("Q is: {:#?}",q);

        let curr_node = q.pop_front().unwrap();

        if curr_node == *end_node_id {
            continue;
        }

        visited_node_id_set.insert(curr_node);
        let current_handle = Handle::pack(curr_node, false);

        //println!("Visited node id set is: {:#?}",visited_node_id_set);
        //println!("Q is: {:#?}",q);
        //println!("Curr_node_is: {}",curr_node);
        //io::stdin().read_line(&mut String::new());

        //Get all paths that end in curr_node
        //let curr_paths_list : Vec<_> = all_paths_list.iter()
        //                                    .filter(|x| x.ends_with(&[curr_node])).collect();
        //let v = Vec::from_iter(curr_paths_list);

        //Get all paths that end in curr_node
        let mut curr_paths_list: Vec<_> = all_paths_list.clone();
        curr_paths_list.retain(|x| x.ends_with(&[curr_node]));

        //Only keep those which don't
        all_paths_list.retain(|x| !x.ends_with(&[curr_node]));

        //println!("Curr_paths_list: {:#?}",curr_paths_list);
        //io::stdin().read_line(&mut String::new());

        for neighbor in handle_edges_iter(g, current_handle, Direction::Right) {
            //Append, for each current_path, this neighbor
            let mut temp = curr_paths_list.clone();
            temp.iter_mut().for_each(|x| x.push(neighbor.id()));
            all_paths_list.append(&mut temp);

            //Add new node to queue
            if !visited_node_id_set.contains(&neighbor.id()) && !q.contains(&neighbor.id()) {
                q.push_back(neighbor.id());
            }
        }

        //println!("All_paths_list: {:#?}",all_paths_list);
        //io::stdin().read_line(&mut String::new());
    }

    //Only keep paths that end in end_node_id
    //start_node_id does not have to be checked
    all_paths_list.retain(|x| x.ends_with(&[*end_node_id]));

    //println!("All paths between {} and {} are: {:#?}",start_node_id, end_node_id, all_paths_list);

    //io::stdin().read_line(&mut String::new());

    all_paths_list
}

/// Detects variants from a list of bubbles
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

        //println!("BEFORE DETECT");

        detect_variants_per_reference(
            &current_ref,
            &ref_path,
            possible_bubbles_list,
            graph,
            node_id_to_path_and_pos_map,
            &mut stuff_to_alts_map,
            verbose,
        );

        //println!("AFTER DETECT");
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
    //println!("BEFORE GET LAST");
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

        //println!("BEFORE FIND START");

        let start_node_index_in_ref_path: usize;
        match ref_path.iter().position(|&r| NodeId::from(r) == *start) {
            None => continue, //ignore, start not found in ref path
            Some(r) => start_node_index_in_ref_path = r,
        };

        //println!("BEFORE PRINT ALL PATHS");

        let all_path_list: Vec<Vec<NodeId>> = find_all_paths_between(&graph, start, end);

        //println!("AFTER PRINT ALL PATHS");

        //println!("All paths list: {:?}",all_path_list);
        for path in &all_path_list {
            if verbose {
                println!("\tPath: {:?}", path);
            }

            //println!("INSIDE FOR LOOP");

            let mut pos_ref = node_id_to_path_and_pos_map[start][current_ref] + 1;
            let mut pos_path = pos_ref;

            let max_index = cmp::min(path.len(), ref_path.len());

            let mut current_index_step_path = 0;
            let mut current_index_step_ref = 0;

            for _i in 0..max_index {
                //Check if ref_path goes out of bounds
                //TODO: check how paths_to_steps is created, there may be some problems there
                // since ref_path is obtained from paths_to_steps
                if current_index_step_ref + start_node_index_in_ref_path >= ref_path.len() {
                    continue;
                }

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

    for (path_name, steps_list) in path_to_steps_map {
        let mut pos = 0;

        for node_id_is_rev in steps_list {
            // Get orientation
            let _is_rev = node_id_is_rev.pop().unwrap();
            // Get the id of the node string -> NodeId
            let node_id: NodeId = NodeId::from(node_id_is_rev.parse::<u64>().unwrap());

            let node_handle = Handle::pack(node_id, false);
            let seq = graph.sequence(node_handle);

            node_id_to_path_and_pos_map
                .entry(node_id)
                .or_insert(HashMap::new());

            if !node_id_to_path_and_pos_map[&node_id].contains_key(path_name) {
                node_id_to_path_and_pos_map
                    .get_mut(&node_id)
                    .unwrap()
                    .insert(String::from(path_name), pos);
            }

            pos += seq.len();
        }
    }

    node_id_to_path_and_pos_map
}

/// Returns the json of a given graph
fn graph_to_json(graph: &HashGraph) -> JsonValue {
    //Create empty json
    let mut graph_json = JsonValue::new_object();

    //Obtains nodes and edges in the graph
    let mut nodes: Vec<u64> = handles_iter(graph).map(|x| u64::from(x.id())).collect();
    let mut edges: Vec<(u64, u64)> = std::iter::from_fn(graph.edges_iter_impl())
        .map(|edge| (u64::from(edge.0.id()), u64::from(edge.1.id())))
        .collect();

    //Sort both vecs so that they're easier to read
    nodes.sort();
    edges.sort();

    //Create an array of objects with from/where for each edge
    let mut edges_array = JsonValue::new_array();
    for (start, end) in edges {
        let temp = object! {
            start: start,
            end: end
        };
        edges_array.push(temp).unwrap();
    }

    graph_json["nodes"] = nodes.into();
    graph_json["edges"] = edges_array;

    graph_json
}

/// Writes the json to a file
fn json_to_file(json: &JsonValue, path: &PathBuf) -> std::io::Result<()> {
    let mut buffer = BufWriter::new(File::create(path)?);
    json.write_pretty(&mut buffer, 4)?;
    Ok(())
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
        .arg(
            Arg::with_name("json")
                .short("j")
                .long("json")
                .value_name("FILE")
                .takes_value(true)
                .help("Sets the path where to store the json of both the starting graph and its bfs-tree"),
        )
        .get_matches();

    let in_path_file = matches
        .value_of("input")
        .expect("Could not parse argument --input");
    let out_path_file = matches
        .value_of("output")
        .expect("Could not parse argument --output");
    let verbose = matches.is_present("verbose");
    let json_out = matches.value_of("json");

    //let in_path_file = "./input/samplePath3.gfa";
    //let out_path_file = "./input/samplePath3.vcf";

    if let Some(gfa) = gfa::parser::parse_gfa(&PathBuf::from(in_path_file)) {
        let graph = HashGraph::from_gfa(&gfa);

        // Stores json of starting graph in a file
        if json_out.is_some() {
            let json = graph_to_json(&graph);
            let buffer = PathBuf::from(format!("{}-graph.json", json_out.unwrap()));
            json_to_file(&json, &buffer).expect("Cannot write json");
        }

        // Obtains, for each path, a list of all its steps (its nodes)
        let mut path_to_steps_map: HashMap<String, Vec<String>> = paths_to_steps(&graph);

        // Obtains, for each node, its position in each path where the node is in
        let node_id_to_path_and_pos_map: BTreeMap<NodeId, HashMap<String, usize>> =
            get_node_positions_in_paths(&graph, &mut path_to_steps_map);

        if verbose {
            for node_id in node_id_to_path_and_pos_map.keys() {
                let path_and_pos_map = node_id_to_path_and_pos_map.get(node_id);
                println!("Node_id : {}", node_id);

                for (path, pos) in path_and_pos_map.unwrap() {
                    println!("Path: {}  -- Pos: {}", path, pos);
                }
            }
        }

        //Obtains the tree representing the bfs
        let g_bfs: HashGraph = bfs(&graph, &NodeId::from(1));

        // Stores json of bfs-tree in a file
        if json_out.is_some() {
            let json = graph_to_json(&g_bfs);
            let buffer = PathBuf::from(format!("{}-bfs.json", json_out.unwrap()));
            json_to_file(&json, &buffer).expect("Cannot write json");
        }

        if verbose {
            for h in handles_iter(&g_bfs) {
                display_node_edges(&g_bfs, &h);
            }
        }

        // Obtains, for each level of the tree, how many nodes are there
        let (distances_map, ordered_node_id_list) = bfs_distances(&g_bfs, &NodeId::from(1));

        if verbose {
            println!("\nNode --> Distance from root");
            for (node_id, distance) in distances_map.iter() {
                println!("{} - distance from root: {}", node_id, distance);
            }
        }

        //Obtains a map where, for each distance from root, the number of nodes
        //at that distance are present
        let dist_to_num_nodes: BTreeMap<u64, usize> = get_dist_to_num_nodes(&distances_map);

        if verbose {
            println!("\nDistance from root --> Num. nodes");
            for (k, v) in dist_to_num_nodes.iter() {
                println!("{} --> {}", k, v);
            }
        }

        let possible_bubbles_list: Vec<(NodeId, NodeId)> =
            detect_bubbles(&distances_map, &ordered_node_id_list, &dist_to_num_nodes);

        if verbose {
            println!("\nBubbles");
            println!("Detected bubbles {:#?}", possible_bubbles_list);
            println!("\n------------------");

            //Obtains, for each path, a string representing all the bases of the path (not actually used)
            let path_to_sequence_map: HashMap<String, String> =
                get_path_to_sequence(&graph, &path_to_steps_map);
            println!("Path to sequence: {:?}", path_to_sequence_map);
        }

        let vcf_list = detect_all_variants(
            &path_to_steps_map,
            &possible_bubbles_list,
            &graph,
            &node_id_to_path_and_pos_map,
            verbose,
        );

        //Write variants to file
        write_to_file(&PathBuf::from(out_path_file), &vcf_list).unwrap();
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

#[cfg(test)]
mod tests {
    use super::*;
    use handlegraph::handlegraph::HandleGraph;
    use handlegraph::hashgraph::HashGraph;

    //Used in other tests
    fn read_test_gfa() -> HashGraph {
        use gfa::parser::parse_gfa;

        HashGraph::from_gfa(&parse_gfa(&PathBuf::from("./input/samplePath3.gfa")).unwrap())
    }

    #[test]
    fn test_path_to_steps() {
        let graph = read_test_gfa();

        let path_to_steps_map: HashMap<String, Vec<String>> = paths_to_steps(&graph);

        // Check if all paths have been found
        assert_eq!(path_to_steps_map.keys().len(), 3);
        assert!(path_to_steps_map.contains_key("x"));
        assert!(path_to_steps_map.contains_key("y"));
        assert!(path_to_steps_map.contains_key("z"));

        //Check for each path that all its node have been found

        //P	x	1+,3+,5+,6+,8+,9+,11+,12+,13+,15+,16+,18+,19+
        let path_x: Vec<String> = vec![
            "1+".to_string(),
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
            "19+".to_string(),
        ];
        assert_eq!(*path_to_steps_map.get("x").unwrap(), path_x);

        //P	y	1+,2+,5+,6+,8+,9+,10+,11+,13+,14+,16+,17+,19+
        let path_y: Vec<String> = vec![
            "1+".to_string(),
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
            "19+".to_string(),
        ];
        assert_eq!(*path_to_steps_map.get("y").unwrap(), path_y);

        //P	z	1+,2+,4+,6+,7+,9+,10+,11+,12+,13+,14+,16+,18+,19+
        let path_z: Vec<String> = vec![
            "1+".to_string(),
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
            "19+".to_string(),
        ];
        assert_eq!(*path_to_steps_map.get("z").unwrap(), path_z);
    }

    #[test]
    fn test_bfs() {
        let graph = read_test_gfa();
        let g_bfs = bfs(&graph, &NodeId::from(1));

        // All nodes must be present in bfs
        assert_eq!(graph.node_count(), g_bfs.node_count());

        // There should be less (or the same number of) edges in g_bfs
        // since nodes can only get added once
        assert!(graph.edge_count() >= g_bfs.edge_count());
    }

    #[test]
    fn test_bfs_distances() {
        let graph = read_test_gfa();
        let g_bfs = bfs(&graph, &NodeId::from(1));

        let (distances_map, ordered_node_id_list) = bfs_distances(&g_bfs, &NodeId::from(1));

        // Should pass easily
        assert_eq!(
            distances_map[&NodeId::from(2)],
            distances_map[&NodeId::from(3)]
        );
        assert_eq!(
            distances_map[&NodeId::from(4)],
            distances_map[&NodeId::from(5)]
        );
        assert_eq!(
            distances_map[&NodeId::from(7)],
            distances_map[&NodeId::from(8)]
        );

        // Caused problems in dfs
        assert_eq!(
            distances_map[&NodeId::from(10)],
            distances_map[&NodeId::from(11)]
        );
        assert_eq!(
            distances_map[&NodeId::from(12)],
            distances_map[&NodeId::from(13)]
        );

        //Should pass easily
        assert_eq!(
            distances_map[&NodeId::from(14)],
            distances_map[&NodeId::from(15)]
        );
        assert_eq!(
            distances_map[&NodeId::from(17)],
            distances_map[&NodeId::from(18)]
        );

        // Check ordered_node_list now
        let mut sorted = ordered_node_id_list.clone();
        sorted.sort();
        assert!(ordered_node_id_list.eq(&sorted));
    }

    #[test]
    fn test_dist_to_num_nodes() {
        let graph = read_test_gfa();
        let g_bfs = bfs(&graph, &NodeId::from(1));

        let (distances_map, _) = bfs_distances(&g_bfs, &NodeId::from(1));
        let mut dist_to_num_nodes: BTreeMap<u64, usize> = BTreeMap::new();
        for (_, distance) in distances_map.iter() {
            if !dist_to_num_nodes.contains_key(&distance) {
                dist_to_num_nodes.insert(*distance, 0);
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
        let g_bfs = bfs(&graph, &NodeId::from(1));

        let (distances_map, ordered_node_id_list) = bfs_distances(&g_bfs, &NodeId::from(1));
        let mut dist_to_num_nodes: BTreeMap<u64, usize> = BTreeMap::new();
        for (_, distance) in distances_map.iter() {
            if !dist_to_num_nodes.contains_key(&distance) {
                dist_to_num_nodes.insert(*distance, 0);
            }
            *dist_to_num_nodes.get_mut(distance).unwrap() += 1;
        }

        let possible_bubbles_list: Vec<(NodeId, NodeId)> =
            detect_bubbles(&distances_map, &ordered_node_id_list, &dist_to_num_nodes);

        //println!("Possible bubbles list: {:#?}",possible_bubbles_list);

        assert!(possible_bubbles_list.contains(&(NodeId::from(1), NodeId::from(6))));
        assert!(possible_bubbles_list.contains(&(NodeId::from(6), NodeId::from(9))));
        assert!(possible_bubbles_list.contains(&(NodeId::from(9), NodeId::from(16))));
        assert!(possible_bubbles_list.contains(&(NodeId::from(16), NodeId::from(19))));
    }

    fn run_whole_script(graph : HashGraph) -> Vec<Variant> {
        //Obtain preliminary data required for future steps
        let mut path_to_steps_map: HashMap<String, Vec<String>> = paths_to_steps(&graph);
        let node_id_to_path_and_pos_map: BTreeMap<NodeId, HashMap<String, usize>> =
            get_node_positions_in_paths(&graph, &mut path_to_steps_map);

        //Compute bfs and analyze the results
        let g_bfs: HashGraph = bfs(&graph, &NodeId::from(1));
        let (distances_map, ordered_node_id_list) = bfs_distances(&g_bfs, &NodeId::from(1));
        let dist_to_num_nodes: BTreeMap<u64, usize> = get_dist_to_num_nodes(&distances_map);

        //Find the bubbles
        let possible_bubbles_list: Vec<(NodeId, NodeId)> =
            detect_bubbles(&distances_map, &ordered_node_id_list, &dist_to_num_nodes);

        //Find variants from bubbles
        let vcf_list = detect_all_variants(
            &path_to_steps_map,
            &possible_bubbles_list,
            &graph,
            &node_id_to_path_and_pos_map,
            false,
        );

        vcf_list
    }

    #[test]
    fn test_variant_detection() {
        let graph = read_test_gfa();
        let variants_found = run_whole_script(graph);

        //Check that all variants have been found
        assert_eq!(variants_found.len(), 21);

        //Check one variant per type

        //x	9	.	G	A	.	.	TYPE=snv	GT	0|1
        let snv = Variant {
            chromosome: "x".to_string(),
            position: "9".to_string(),
            id: ".".to_string(),
            reference: "G".to_string(),
            alternate: "A".to_string(),
            quality: ".".to_string(),
            filter: ".".to_string(),
            info: format!("TYPE={}", "snv"),
            format: "GT".to_string(),
            sample_name: "0|1".to_string(),
        };
        assert!(variants_found.contains(&snv));

        //x	18	.	T	TAA	.	.	TYPE=ins	GT	0|1
        let ins = Variant {
            chromosome: "x".to_string(),
            position: "18".to_string(),
            id: ".".to_string(),
            reference: "T".to_string(),
            alternate: "TAA".to_string(),
            quality: ".".to_string(),
            filter: ".".to_string(),
            info: format!("TYPE={}", "ins"),
            format: "GT".to_string(),
            sample_name: "0|1".to_string(),
        };
        assert!(variants_found.contains(&ins));

        //y	18	.	TAA	T	.	.	TYPE=del	GT	0|1
        let del = Variant {
            chromosome: "y".to_string(),
            position: "18".to_string(),
            id: ".".to_string(),
            reference: "TAA".to_string(),
            alternate: "T".to_string(),
            quality: ".".to_string(),
            filter: ".".to_string(),
            info: format!("TYPE={}", "del"),
            format: "GT".to_string(),
            sample_name: "0|1".to_string(),
        };
        assert!(variants_found.contains(&del));
    }

    #[test]
    fn find_all_paths_1() {
       let mut graph = HashGraph::new();
       
       //Add nodes
       let h1 = graph.append_handle("A");
       let h2 = graph.append_handle("T");
       let h3 = graph.append_handle("C");
       let h4 = graph.append_handle("G");
       let h5 = graph.append_handle("AC");

       //Add edges
       graph.create_edge(&Edge(h1,h2));
       graph.create_edge(&Edge(h2,h3));
       graph.create_edge(&Edge(h3,h4));
       //Loop
       graph.create_edge(&Edge(h3,h5));
       graph.create_edge(&Edge(h5,h3));

       let paths = find_all_paths_between(&graph, &h1.id(), &h4.id());

       assert!(paths.len()==1);
       assert!(paths.contains(&vec![h1.id(),h2.id(),h3.id(),h4.id()]));
    }

    #[test]
    fn find_all_paths_2() {
       let mut graph = HashGraph::new();
       
       //Add nodes
       let h1 = graph.append_handle("A");
       let h2 = graph.append_handle("T");
       let h3 = graph.append_handle("C");
       let h4 = graph.append_handle("G");

       //Add edges
       //Path 1
       graph.create_edge(&Edge(h1,h2));
       //Path 2
       graph.create_edge(&Edge(h1,h3));
       graph.create_edge(&Edge(h3,h2));
       //Path 3
       graph.create_edge(&Edge(h3,h4));
       graph.create_edge(&Edge(h4,h2));
    
       let paths = find_all_paths_between(&graph, &h1.id(), &h2.id());

       assert!(paths.len()==3);
       assert!(paths.contains(&vec![h1.id(),h2.id()]));
       assert!(paths.contains(&vec![h1.id(),h3.id(),h2.id()]));
       assert!(paths.contains(&vec![h1.id(),h3.id(),h4.id(),h2.id()]));
    }

    #[test]
    fn find_bfs_1() {
       let mut graph = HashGraph::new();
       
       //Add nodes
       let h1 = graph.append_handle("A");
       let h2 = graph.append_handle("T");
       let h3 = graph.append_handle("C");
       let h4 = graph.append_handle("G");
       let h5 = graph.append_handle("AC");

       //Add edges
       graph.create_edge(&Edge(h1,h2));
       graph.create_edge(&Edge(h2,h3));
       graph.create_edge(&Edge(h3,h4));
       //Loop
       graph.create_edge(&Edge(h3,h5));
       graph.create_edge(&Edge(h5,h3));

       let g_bfs = bfs(&graph, &NodeId::from(1));

       //Check g_bfs does not contain the loop
       assert!(g_bfs.has_edge(h3,h5));
       assert!(!g_bfs.has_edge(h5,h3));
    }

    #[test]
    fn find_bfs_paths_2() {
       let mut graph = HashGraph::new();
       
       //Add nodes
       let h1 = graph.append_handle("A");
       let h2 = graph.append_handle("T");
       let h3 = graph.append_handle("C");
       let h4 = graph.append_handle("G");

       //Add edges
       //Path 1
       graph.create_edge(&Edge(h1,h2));
       //Path 2
       graph.create_edge(&Edge(h1,h3));
       graph.create_edge(&Edge(h3,h2));
       //Path 3
       graph.create_edge(&Edge(h3,h4));
       graph.create_edge(&Edge(h4,h2));
    
       let g_bfs = bfs(&graph, &NodeId::from(1));

        assert!(g_bfs.has_edge(h1, h2));
        assert!(g_bfs.has_edge(h1,h3));
        assert!(g_bfs.has_edge(h3,h4));

        //These will represent bubbles
        assert!(!g_bfs.has_edge(h3, h2));
        assert!(!g_bfs.has_edge(h4, h2));      
    }

    #[test]
    fn find_bubbles_1() {
       let mut graph = HashGraph::new();
       
       //Add nodes
       let h1 = graph.append_handle("A");
       let h2 = graph.append_handle("T");
       let h3 = graph.append_handle("C");
       let h4 = graph.append_handle("G");
       let h5 = graph.append_handle("AC");

        //Add edges
        graph.create_edge(&Edge(h1,h2));
        graph.create_edge(&Edge(h2,h3));
        graph.create_edge(&Edge(h3,h4));
        //Loop
        graph.create_edge(&Edge(h3,h5));
        graph.create_edge(&Edge(h5,h3));

       let g_bfs = bfs(&graph, &NodeId::from(1));

       let (distances_map, ordered_node_id_list) = bfs_distances(&g_bfs, &NodeId::from(1));
       let mut dist_to_num_nodes: BTreeMap<u64, usize> = BTreeMap::new();
       for (_, distance) in distances_map.iter() {
           if !dist_to_num_nodes.contains_key(&distance) {
               dist_to_num_nodes.insert(*distance, 0);
           }
           *dist_to_num_nodes.get_mut(distance).unwrap() += 1;
       }

       let possible_bubbles_list: Vec<(NodeId, NodeId)> =
           detect_bubbles(&distances_map, &ordered_node_id_list, &dist_to_num_nodes);

       //println!("Possible bubbles list: {:#?}",possible_bubbles_list);   

       assert!(possible_bubbles_list.is_empty());
    }

    #[test]
    fn find_bubbles_2() {
        let mut graph = HashGraph::new();
       
        //Add nodes
        let h1 = graph.append_handle("A");
        let h2 = graph.append_handle("T");
        let h3 = graph.append_handle("C");
        let h4 = graph.append_handle("G");

        //Add edges
        //Path 1
        graph.create_edge(&Edge(h1,h2));
        //Path 2
        graph.create_edge(&Edge(h1,h3));
        graph.create_edge(&Edge(h3,h2));
        //Path 3
        graph.create_edge(&Edge(h3,h4));
        graph.create_edge(&Edge(h4,h2));
    
        let g_bfs = bfs(&graph, &NodeId::from(1));

        let (distances_map, ordered_node_id_list) = bfs_distances(&g_bfs, &NodeId::from(1));
        let mut dist_to_num_nodes: BTreeMap<u64, usize> = BTreeMap::new();
        for (_, distance) in distances_map.iter() {
            if !dist_to_num_nodes.contains_key(&distance) {
                dist_to_num_nodes.insert(*distance, 0);
            }
            *dist_to_num_nodes.get_mut(distance).unwrap() += 1;
        }
 
        let possible_bubbles_list: Vec<(NodeId, NodeId)> =
            detect_bubbles(&distances_map, &ordered_node_id_list, &dist_to_num_nodes);
 
        //println!("Possible bubbles list: {:#?}",possible_bubbles_list);   
 
        assert!(possible_bubbles_list.contains(&(NodeId::from(1), NodeId::from(2)))); 
    }


}
