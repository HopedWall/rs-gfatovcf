use std::collections::HashMap;
//use std::collections::BTreeMap; //like hashmap but sorted
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

fn process_step(g : &HashGraph, s : &PathStep) -> String {
    let h = g.get_handle_of_step(s);
    let is_rev = h.is_reverse();
    let id = h.id();
    
    //Added display trait in rs-handlegraph
    let result = id.to_string();
    
    if is_rev == true {
        result.push_str("+");
    } else {
        result.push_str("-");
    }

    result
}

fn create_into_hashmap(g : &HashGraph, path_to_steps : &HashMap<String, Vec<String>>,  path : &PathId, step : &PathStep) {
    let path_name = g.get_path_name(path);
    if !path_to_steps.contains_key(path_name) {
        path_to_steps[path_name] = Vec::new();
    } 
    path_to_steps[path_name].push(process_step(g,step));
}


// TODO: fix last argument
fn create_edge_and_so_on(g : &HashGraph, g_dfs : &HashGraph, handle1 : &Handle, handle2 : &Handle, so_on_function : &dyn Fn(&HashGraph, &Handle), args : &Vec<String>) {
    let handle1_id = g.get_id(handle1);
    let handle2_id = g.get_id(handle2); 

    if !g_dfs.has_node(handle2_id) {
        so_on_function(args[0],args[1]);
        
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
//TODO: see create_edge_and_so_on
fn dfs(g : &HashGraph, g_dfs : &HashGraph, node_id : NodeId) {
    let current_node = g.get_handle(node_id, false);
    let sequence_node = g.get_sequence(&current_node);

    //TODO: fix this
    g.follow_edges(
        &current_node, 
        Direction::Right,  //What should go here?
        |neighbor| {
            create_edge_and_so_on(
                &g, &g_dfs, &current_node, &neighbor, &dfs(g, g_dfs, neighbor), g.get_id(neighbor));
            true
        });
}

fn calculate_distance(visited_node_id_set : HashSet<String>, prev_node_id : String, neighbour_id : String, Q : VecDeque<String>, distances_map : HashMap<String,u64>) {
    if !visited_node_id_set.contains(&neighbour_id) {
        distances_map[&neighbour_id] = distances_map[&prev_node_id] + 1;
        // is push_back correct?
        Q.push_back(neighbour_id);
        visited_node_id_set.insert(neighbour_id);
    }
}

// This function can be optimized via closures!
fn bfs_distances(g : &HashGraph, starting_node_id : &PathId, distances_map : &HashMap<u64, u64>) {
    let mut visited_node_id_set : HashSet<u64> = HashSet::new();
    let node_id_list = Vec::new();
    g.for_each_handle(|h| true);

    for node_id in node_id_list {
        distances_map[node_id] = 0;
    }
    node_id_list.clear();

    //TODO: complete this function
}

fn show_edge(g_dfs : &HashGraph, a : &Handle, b : &Handle) {
     print!("{} --> {}", g_dfs.get_id(a), g_dfs.get_id(b));
}

fn display_node_edges(g_dfs : &HashGraph, h : &PathId) {
    print!("node {}", g_dfs.get_id(h));
    g_dfs.follow_edges(
        h, Direction::Right,
        |n| {show_edge(&g_dfs,h,n)}
    );
} //TODO: complete this function

//TODO: complete these function
//fn print_all_paths_util()
//fn print_all_paths()

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

            for (path, pos) in path_and_pos_map.iter() {
                println!("Path: {}  -- Pos: {}", path, pos);
            }
        }

        let start_node = graph.get_handle(NodeId::from(1), false);

        let g_dfs = HashGraph::new();

        dfs(&graph, &g_dfs, NodeId::from(1));

        g_dfs.for_each_handle(display_node_edges);

        let value = bfs_distances(&g_dfs, &1);

        let distances_map = value[0];
        let ordered_node_id_list = value[1];

        for (node_id, distance) in distances_map.iter() {
            println!("{} - distance from root: {}", node_id, distance);
        }

        let dist_to_num_nodes = HashMap::new();

        for (node_id, distance) in distances_map.iter() {
            if !distances_map.contains_key(distance) {
                dist_to_num_nodes[distance] = 0;
            }
            dist_to_num_nodes[distance] += 1;
        }

        print!("Distance from root --> Num. nodes");
        for (k,v) in dist_to_num_nodes.iter() {
            println!("{} --> {}", k, v);
        }

        println!("Bubbles");
        let possible_bubbles_list = Vec::new();
        let first_bubble = true;

        for node_id in ordered_node_id_list {
            let key = distances_map[node_id];
            if dist_to_num_nodes[key] == 1 {
                if !first_bubble {
                    //println!("{} END {} {}", node_id, node_id_to_path_and_pos_map[node_id],g_dfs.get_sequence(g_dfs.get_handle(node_id)));
                    // REMOVE value here
                }
                first_bubble = false;
                //println!("{} START {}");
            }     
        } //TODO: complete this

        println!("\n------------------");

        let path_to_sequence_map = HashMap::new();
        for (path_name, steps_list) in path_to_steps_map.iter() {
            path_to_sequence_map[path_name] = String::new();

            for node_id_rev in steps_list {
                path_to_sequence_map[path_name].push_str(graph.get_sequence(graph.get_handle(&node_id_rev, false)));
            }
        }

        let stuff_to_alts_dict = HashMap::new();
        for current_ref in path_to_steps_map.keys() {
            
            let path = Vec::new(); //TODO: fix this variable

            for (start,end) in possible_bubbles_list[0:possible_bubbles_list.len()-1] {
                println!("ref_path: {}",ref_path);
                println!("Bubble [{},{}]",start, end);
                start_node_index_in_ref_path = ref_path.index(start);
                let all_path_list = Vec::new();
                //print_all_paths(g, start, end, all_path_list)

                for path in all_path_list {
                    println!("\tPath: {}", path);
                    let pos_ref = node_id_to_path_and_pos_map[tart][current_ref]+1;
                    let pos_path = pos_ref;

                    println!("Start paths position: {}",pos_ref);

                    let max_index = min(path.len(), ref_path.len());
                    let current_index_step_path = 0;
                    let current_index_step_ref = 0;

                    for i in range(0, max_index) {

                        let current_node_id_path = path[current_index_step_path];
                        let current_node_id_ref = ref_path[current_index_step_ref + start_node_index_in_ref_path];
                    
                        println!("{} {} ---> {} {}", pos_ref, pos_path, current_node_id_ref, current_node_id_path);
                        
                        if current_node_id_ref == current_node_id_path {
                            println!("REFERENCE");
                            let node_seq = g.get_sequence(g.get_handle(current_node_id_ref));
                            pos_ref += len(node_seq);
                            pos_path = pos_ref;
                            current_index_step_ref += 1;
                            current_index_step_path += 1;
                        } else {
                            let succ_node_id_path = path[current_index_step_path + 1];
                            let succ_node_id_ref = ref_path[current_index_step_ref + start_node_index_in_ref_path + 1];
                            if succ_node_id_ref == current_node_id_path {
                                println!("DEL");
                                let node_seq_ref = g.get_sequence(g.get_handle(current_node_id_ref));
                                let prec_node_id_ref = ref_path[current_index_step_ref + start_node_index_in_ref_path - 1];
                                let prec_nod_seq_ref = g.get_sequence(g.get_handle(prec_node_id_ref))
                                //key = '_'.join([current_ref, str(pos_path - 1), prec_nod_seq_ref[-1] + node_seq_ref])
                                if key not in stuff_to_alts_dict {
                                    stuff_to_alts_dict[key] = HashSet::new();
                                }
                                //stuff_to_alts_dict[key].add(prec_nod_seq_ref[-1] + '_del')
                                pos_ref += len(node_seq_ref);
                                current_index_step_ref += 1;
                                current_node_id_ref = ref_path[current_index_step_ref + start_node_index_in_ref_path -1];
                                println!("\t {}", current_node_id_ref);
                                continue;
                            } else if (succ_node_id_path == current_node_id_ref) {
                                println!("INS");
                                let node_seq_path = g.get_sequence(g.get_handle(current_node_id_path))
                                
                                let prec_node_id_ref = ref_path[current_index_step_ref + start_node_index_in_ref_path-1];
                                let prec_nod_seq_ref = g.get_sequence(g.get_handle(prec_node_id_ref));
                                //key = '_'.join([current_ref, str(pos_ref-1), prec_nod_seq_ref[-1]])
                                
                                if key not in stuff_to_alts_dict {
                                    stuff_to_alts_dict[key] = HashSet::new();
                                }
                                //stuff_to_alts_dict[key].add(prec_nod_seq_ref[-1] + node_seq_path + '_ins')
                                pos_path += len(node_seq_path)
                                current_index_step_path += 1
                                current_node_id_path = path[current_index_step_path]
                                println!("\t{}", current_node_id_path)
                                continue;
                            } else {
                                node_seq_ref = g.get_sequence(g.get_handle(current_node_id_ref));
                                node_seq_path = g.get_sequence(g.get_handle(current_node_id_path));

                                if node_seq_ref == node_seq_path {
                                    println!("REFERENCE");
                                } else {
                                    println!("SNV");
                                }

                                key = '_'.join([current_ref, str(pos_path), node_seq_ref])
                                if key not in stuff_to_alts_dict{
                                    stuff_to_alts_dict[key] = set();
                                }
                                    
                                stuff_to_alts_dict[key].add(node_seq_path + "snv");

                                pos_ref += len(node_seq_ref);
                                pos_path += len(node_seq_path);
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
