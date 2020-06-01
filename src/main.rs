use std::collections::HashMap;
//use std::collections::BTreeMap; //like hashmap but sorted
use std::path::PathBuf;
use handlegraph::graph::HashGraph;
use handlegraph::handlegraph::HandleGraph;
use handlegraph::handle::{NodeId,Handle, Direction};
use handlegraph::pathgraph::PathHandleGraph;
use gfa::gfa::GFA;
use std::fs::File;
use std::io::Write;
extern crate clap;
//use clap::{Arg, App};
extern crate chrono;
use chrono::Utc;
use std::io::{BufReader,BufRead};
use petgraph::Graph;
use std::collections::VecDeque;
use std::collections::HashSet;


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

// Why doesn't this work?
fn process_step(g : &HashGraph, s : &HashGraph::StepHandle) -> String {
    let h = g.get_handle_of_step(s);
    let is_rev = h.is_reverse();
    let id = h.id();
    
    // Should be converted to u64 and then String?
    let result = String::from(id);
    
    if is_rev == true {
        result.push_str("+");
    } else {
        result.push_str("-");
    }

    result
}

fn create_into_hashmap(g : &HashGraph, path_to_steps : &HashMap<String, Vec<String>>,  path : &HashGraph::PathHandle, step : &HashGraph::StepHandle) {
    let path_name = g.get_path_name(path);
    if !path_to_steps.contains_key(path_name) {
        path_to_steps[path_name] = Vec::new();
    } 
    path_to_steps[path_name].push(process_step(g,step));
}

// ------------ FROM NOW ON -- ODGI REQUIRED? ---------

// What's the fourth parameter for?
// fn create_edge_and_so_on(g : HashGraph, g_dfs : Vec<String>, handle1 : HashGraph::NodeHandle, handle2 : HashGraph::NodeHandle, so_on_function : &dyn Fn() -> i32, args : Fn()) {
//     let handle1_id = g.get_id(handle1);
//     let handle2_id = g.get_id(handle2); 

//     if !g_dfs.contains(handle1) {
//         g_dfs.append(g.get_sequence(handle1));
//     }

//     if !g_dfs.contains(handle2) {
//         g_dfs.append(g.get_sequence(&handle2));
//     }
// } //TODO: understand this function better

// fn dfs(g : HashGraph, node_id : HashGraph::NodeHandle) {
//     let current_node = g.get_handle(node_id, false);
//     let sequence_node = g.get_sequence(&current_node);

//     g.follow_edges(
//         &current_node, 
//         Direction::Right,  //What should go here?
//         |neighbor| {
//             //create_edge_and_so_on(g, g_dfs: Vec<String>, current_node, neighbor, dfs, g.get_id(neighbor));
//             true
//         });
// }

// fn calculate_distance(visited_node_id_set : HashSet<String>, prev_node_id : String, neighbour_id : String, Q : VecDeque<String>, distances_map : HashMap<String,u64>) {
//     if !visited_node_id_set.contains(&neighbour_id) {
//         distances_map[&neighbour_id] = distances_map[&prev_node_id] + 1;
//         // is push_back correct?
//         Q.push_back(neighbour_id);
//         visited_node_id_set.insert(neighbour_id);
//     }
// }

// // This function can be optimized via closures!
// fn bfs_distances(g : HashGraph, starting_node_id : HashGraph::NodeHandle) {
//     let visited_node_id_set = HashSet::new();
//     let node_id_list = Vec::new();
//     g.for_each_handle(|h| true);

//     for node_id in node_id_list {
//         distances_map[node_id] = 0;
//     }
//     node_id_list.clear();

//     //TODO: complete this function
// }

// // fn show_edge(g_dfs : Graph, a : u64, b : u64) {
    
// //     print!(g_dfs.get_id)
// // }




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

    let ref_path = "./input/samplePath3.fa";
    let in_path = "./input/samplePath3.gfa";
    let out_path = "./input/samplePath3.vcf";
      
    if let Some(gfa) = gfa::parser::parse_gfa(&PathBuf::from(in_path)) {
        
        let graph = HashGraph::from_gfa(&gfa);
        println!("{:?}",graph);

        let path_to_steps_map = HashMap::new();

        let get_path_to_steps = |p| {
            graph.for_each_step_in_path(p, |s| create_into_hashmap(&graph, &path_to_steps_map, &p, &s))
        };

        graph.for_each_handle(get_path_to_steps);    
        
        let node_id_to_path_and_pos_map = HashMap::new();
        for (path_name, steps_list) in path_to_steps_map {
            let pos = 0;

            for nodeId_isRev in steps_list {
                
                // TODO: check what these are
                let node_id = nodeId_isRev.parse::<u64>().unwrap();
                let is_rev = nodeId_isRev;

                // Why does get handle require 2 parameters?
                let node_handle = graph.get_handle(NodeId::from(node_id), false);
                let seq = graph.get_sequence(&node_handle);

                if !node_id_to_path_and_pos_map.contains_key(&node_id) {
                    node_id_to_path_and_pos_map[&node_id] = HashMap::new();
                }

                if !node_id_to_path_and_pos_map[&node_id].contains_key(&path_name) {
                    node_id_to_path_and_pos_map[&node_id][&path_name] = pos;
                }

                pos += seq.len();
            }
        }

        // TODO: this must be sorted
        //sorted = new BTreeMap(node_id_to_path_and_pos_map);
        for node_id in node_id_to_path_and_pos_map.keys() {
            let path_and_pos_map = node_id_to_path_and_pos_map[node_id];
            println!("Node_id : {}", node_id);

            for (path, pos) in path_and_pos_map {
                println!("Path: {}  -- Pos: {}", path, pos);
            }
        }

        let start_node = graph.get_handle(NodeId::from(1), false);

        //let g_dfs = Vec::new(); //Will be used a stack for DFS
        //Not an array as comment states, but a graph/tree
        //let g_dfs = petgraph::Graph::new();

        //dfs(graph,1);

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
