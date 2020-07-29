//! # GFAtoVCF
//! `GFAtoVCF` is a tool that finds variants in a Variation Graph.

use handlegraph::handle::NodeId;
use handlegraph::handlegraph::handles_iter;
use handlegraph::hashgraph::HashGraph;
use std::collections::BTreeMap; //like hashmap but sorted
use std::collections::HashMap;
use std::path::PathBuf;
extern crate clap;
use clap::{App, Arg};

// Import bubble detection functions
mod bubble_detection;
use crate::bubble_detection::*;

// Import variant identification functions
mod variant_identification;
use crate::variant_identification::*;

//Import fileIO
mod file_io;
use crate::file_io::*;

/// The function that runs the script
fn main() {
    let matches = App::new("rs-GFAtoVCF")
        .version("2.0")
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
                .help("Shows debug messages during execution"),
        )
        .arg(
            Arg::with_name("json")
                .short("j")
                .long("json")
                .value_name("PATH")
                .takes_value(true)
                .help("Sets the path where to store the json of both the starting graph and its bfs-tree"),
        )
        .arg(
            Arg::with_name("max-edges")
                .short("m")
                .long("max-edges")
                .value_name("NUMBER")
                .takes_value(true)
                .default_value("100")
                .help("Sets the maximum amount of edges to be used to find paths between nodes"),
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
    let max_edges = matches
        .value_of("max-edges")
        .unwrap()
        .parse::<i32>()
        .unwrap();

    if let Some(gfa) = gfa::parser::parse_gfa(&PathBuf::from(in_path_file)) {
        let graph = HashGraph::from_gfa(&gfa);

        // Stores json of starting graph in a file
        if json_out.is_some() {
            let json = graph_to_json(&graph);
            let buffer = PathBuf::from(format!("{}graph.json", json_out.unwrap()));
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

        // Obtains the tree representing the bfs
        let g_bfs: HashGraph = bfs(&graph, &NodeId::from(1));

        // Stores json of bfs-tree in a file
        if json_out.is_some() {
            let json = graph_to_json(&g_bfs);
            let buffer = PathBuf::from(format!("{}bfs.json", json_out.unwrap()));
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

        // Obtains a map where, for each distance from root, the number of nodes
        // at that distance are present
        let dist_to_num_nodes: BTreeMap<u64, usize> = get_dist_to_num_nodes(&distances_map);

        if verbose {
            println!("\nDistance from root --> Num. nodes");
            for (k, v) in dist_to_num_nodes.iter() {
                println!("{} --> {}", k, v);
            }
        }

        // Obtains a list of bubbles, each Bubble is represented as (Start_Node_Id, End_Node_Id)
        let possible_bubbles_list: Vec<(NodeId, NodeId)> =
            detect_bubbles(&distances_map, &ordered_node_id_list, &dist_to_num_nodes);

        if verbose {
            println!("\nBubbles");
            println!("Detected bubbles {:#?}", possible_bubbles_list);
            println!("\n------------------");

            // Obtains, for each path, a string representing all the bases of the path
            let path_to_sequence_map: HashMap<String, String> =
                get_path_to_sequence(&graph, &path_to_steps_map);
            println!("Path to sequence: {:?}", path_to_sequence_map);
        }

        // Obtains the list of Variants in the graph
        let vcf_list: Vec<Variant> = detect_all_variants(
            &path_to_steps_map,
            &possible_bubbles_list,
            &graph,
            &node_id_to_path_and_pos_map,
            verbose,
            max_edges,
        );

        // Get path names
        let mut paths_list: Vec<String> = graph
            .paths
            .iter()
            .map(|(_, value)| value.name.clone())
            .collect();
        paths_list.sort();
        // Write variants to file
        write_variants_to_file(&PathBuf::from(out_path_file), &mut paths_list, &vcf_list).unwrap();
    } else {
        panic!("Couldn't parse gfa file!");
    }
}
