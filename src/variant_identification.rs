use handlegraph::handle::{Direction, Handle, NodeId};
use handlegraph::handlegraph::handle_edges_iter;
use handlegraph::handlegraph::HandleGraph;
use handlegraph::hashgraph::HashGraph;
use std::collections::BTreeMap; //like hashmap but sorted
use std::collections::HashMap;
extern crate chrono;
use crate::bubble_detection::Bubble;
use log::info;
use std::cmp;
use std::cmp::Ordering;
use std::collections::HashSet;
use std::collections::VecDeque;

/// A struct that holds Variants, as defined in the VCF format
#[derive(Debug, PartialEq)]
pub struct Variant {
    chromosome: String,
    position: i32,
    id: Option<String>,
    reference: String,
    alternate: Option<String>,
    quality: Option<i32>,
    filter: Option<String>,
    info: Option<String>,
    format: Option<String>,
    sample_name: Option<String>,
}

impl Variant {
    // Returns a Variant as a line of a VCF file
    pub fn to_vcf_string(&self) -> String {
        // Extract quality from Variant
        // Note: quality is of type Option<i32>, when it is None "." should be returned
        let quality_string: String;
        if self.quality.is_none() {
            quality_string = ".".to_string();
        } else {
            quality_string = self.quality.unwrap().to_string();
        }

        let vcf_string = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            self.chromosome,
            self.position,
            self.id.clone().unwrap_or(".".to_string()),
            self.reference,
            self.alternate.clone().unwrap_or(".".to_string()),
            quality_string,
            self.filter.clone().unwrap_or(".".to_string()),
            self.info.clone().unwrap_or(".".to_string()),
            self.format.clone().unwrap_or(".".to_string()),
            self.sample_name.clone().unwrap_or(".".to_string())
        );

        vcf_string
    }
}

/// Detects variants from a list of bubbles
pub fn detect_all_variants(
    path_to_steps_map: &HashMap<String, Vec<String>>,
    possible_bubbles_list: &[Bubble],
    graph: &HashGraph,
    node_id_to_path_and_pos_map: &BTreeMap<NodeId, HashMap<String, usize>>,
    verbose: bool,
    max_edges: i32,
    reference_paths: &Vec<String>,
) -> Vec<Variant> {
    let mut stuff_to_alts_map: HashMap<String, HashSet<String>> = HashMap::new();

    // For each reference path, explore all bubbles in order to find variants;
    // these will be stored in stuff_to_alts_map
    for current_ref in reference_paths {
        // Obtain all steps for current_ref
        let ref_path: Vec<u64> = path_to_steps_map[current_ref]
            .iter()
            .map(|x| x.parse::<u64>().unwrap())
            .collect();

        if verbose {
            println!("path_to_steps_map: {:?}", path_to_steps_map);
        }

        info!("BEFORE DETECT");

        // Loop through all bubbles in order to find all variants
        // for the given reference
        detect_variants_per_reference(
            &current_ref,
            &ref_path,
            possible_bubbles_list,
            graph,
            node_id_to_path_and_pos_map,
            &mut stuff_to_alts_map,
            verbose,
            max_edges,
        );

        info!("AFTER DETECT");
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
        let mut types = "TYPE=".to_string();
        types.push_str(&type_set.join(";TYPE="));

        let v = Variant {
            chromosome: chrom.to_string(),
            position: pos.to_string().parse::<i32>().unwrap(),
            id: None,
            reference: refr.to_string(),
            alternate: Some(alts),
            quality: None,
            filter: None,
            info: Some(types),
            format: Some("GT".to_string()),
            sample_name: Some("0|1".to_string()),
        };

        vcf_list.push(v);
    }

    // Sort vcf_list for printing variants in the correct order
    vcf_list.sort_by(|a, b| match a.chromosome.cmp(&b.chromosome) {
        Ordering::Equal => a.position.cmp(&b.position),
        other => other,
    });

    vcf_list
}
/// Detect variants for a specific reference
fn detect_variants_per_reference(
    current_ref: &str,
    ref_path: &[u64],
    possible_bubbles_list: &[Bubble],
    graph: &HashGraph,
    node_id_to_path_and_pos_map: &BTreeMap<NodeId, HashMap<String, usize>>,
    stuff_to_alts_map: &mut HashMap<String, HashSet<String>>,
    verbose: bool,
    max_edges: i32,
) {
    info!("BEFORE GET LAST");

    // Create closure that will be used later
    let get_last = |prec_node_seq_ref: &str, node_seq_ref| {
        let mut last = prec_node_seq_ref[prec_node_seq_ref.len() - 1..].to_string();
        last.push_str(node_seq_ref);
        last
    };

    // Check all bubbles
    for bubble in possible_bubbles_list {
        let start = bubble.start;
        let end = bubble.end;

        if verbose {
            println!("ref_path: {:?}", ref_path);
            println!("Bubble [{},{}]", start, end);
        }

        info!("BEFORE FIND START");

        let start_node_index_in_ref_path: usize;
        match ref_path.iter().position(|&r| NodeId::from(r) == start) {
            None => continue, //ignore, start not found in ref path
            Some(r) => start_node_index_in_ref_path = r,
        };

        info!("BEFORE FIND ALL PATHS BETWEEN");

        let all_path_list: Vec<Vec<NodeId>> =
            find_all_paths_between(&graph, &start, &end, max_edges);

        info!("AFTER FIND ALL PATHS BETWEEN");

        info!("All paths list: {:?}", all_path_list);
        for path in &all_path_list {
            if verbose {
                println!("\tPath: {:?}", path);
            }

            //println!("INSIDE FOR LOOP");

            let mut pos_ref = node_id_to_path_and_pos_map[&start][current_ref] + 1;
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

/// Return a list of paths between two given nodes, where each path is represented
/// as a list of NodeIds
pub fn find_all_paths_between(
    g: &HashGraph,
    start_node_id: &NodeId,
    end_node_id: &NodeId,
    max_edges: i32,
) -> Vec<Vec<NodeId>> {
    let mut all_paths_list: Vec<Vec<NodeId>> = Vec::new();

    // Put a limit on the maximum amount of edges that can be traversed
    // this should prevent eccessive memory usage
    info!("Max edges is {:#?}", max_edges);
    let mut curr_edges = 0;
    let mut edges_limit_reached = false;

    // Keep a set of visited nodes so that loops are avoided
    let mut visited_node_id_set: HashSet<NodeId> = HashSet::new();

    // Create queue
    // NOTE: this is a Queue based implementation, this was done
    // in order not to get a stack overflow (the previous recursion-based
    // version was often experiencing this kind of issue)
    let mut q: VecDeque<NodeId> = VecDeque::new();

    // Insert first value
    q.push_back(*start_node_id);
    all_paths_list.push(vec![*start_node_id]);

    while !q.is_empty() {
        info!("All paths is {:#?}", all_paths_list);
        info!("Q is: {:#?}", q);

        let curr_node = q.pop_front().unwrap();
        info!("Curr node is {:#?}", curr_node);

        if curr_node == *end_node_id {
            continue;
        }

        visited_node_id_set.insert(curr_node);
        let current_handle = Handle::pack(curr_node, false);

        // Get all paths that end in curr_node
        let mut curr_paths_list: Vec<_> = all_paths_list.clone();
        curr_paths_list.retain(|x| x.ends_with(&[curr_node]));

        // Only keep those which don't
        all_paths_list.retain(|x| !x.ends_with(&[curr_node]));

        info!("Curr_paths_list: {:#?}", curr_paths_list);
        //io::stdin().read_line(&mut String::new());

        for neighbor in handle_edges_iter(g, current_handle, Direction::Right) {
            info!("Neighbor: {:#?}",neighbor.id());
            // Append, for each current_path, this neighbor
            let mut temp = curr_paths_list.clone();
            temp.iter_mut().for_each(|x| x.push(neighbor.id()));
            all_paths_list.append(&mut temp);

            // Add new node to queue
            if !visited_node_id_set.contains(&neighbor.id()) && !q.contains(&neighbor.id()) {
                q.push_back(neighbor.id());
            }

            // Break if too many edges have been visited
            curr_edges = curr_edges + 1;
            if curr_edges > max_edges {
                edges_limit_reached = true;
                break;
            }
        }

        if edges_limit_reached {
            break;
        }

        info!("All_paths_list: {:#?}", all_paths_list);
        //io::stdin().read_line(&mut String::new());
    }

    // Only keep paths that end in end_node_id
    // start_node_id does not have to be checked
    // TODO: maybe not needed?
    all_paths_list.retain(|x| x.ends_with(&[*end_node_id]));

    info!(
        "All paths between {} and {} are: {:#?}",
        start_node_id, end_node_id, all_paths_list
    );

    //io::stdin().read_line(&mut String::new());

    all_paths_list
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bubble_detection::*;
    use handlegraph::handle::{Edge, NodeId};
    use handlegraph::hashgraph::HashGraph;
    use handlegraph::mutablehandlegraph::*;
    use std::path::PathBuf;
    use handlegraph::handlegraph::handles_iter;
    use handlegraph::hashgraph::Path;
    use handlegraph::pathgraph::PathHandleGraph;

    //Used in other tests
    fn read_test_gfa() -> HashGraph {
        use gfa::parser::parse_gfa;

        HashGraph::from_gfa(&parse_gfa(&PathBuf::from("./input/samplePath3.gfa")).unwrap())
    }

    fn run_whole_script(graph: HashGraph) -> (Vec<Variant>, HashGraph, Vec<Bubble>) {
        //Obtain preliminary data required for future steps
        let mut path_to_steps_map: HashMap<String, Vec<String>> = paths_to_steps(&graph);
        let node_id_to_path_and_pos_map: BTreeMap<NodeId, HashMap<String, usize>> =
            get_node_positions_in_paths(&graph, &mut path_to_steps_map);

        //Compute bfs and analyze the results
        let g_bfs: HashGraph = bfs(&graph, &NodeId::from(1));
        let (distances_map, ordered_node_id_list) = bfs_distances(&g_bfs, &NodeId::from(1));
        let dist_to_num_nodes: BTreeMap<u64, usize> = get_dist_to_num_nodes(&distances_map);

        //Find the bubbles
        let possible_bubbles_list: Vec<Bubble> =
            detect_bubbles(&distances_map, &ordered_node_id_list, &dist_to_num_nodes);

        let paths_list = path_to_steps_map
            .keys()
            .map(|x| x.to_string())
            .collect::<Vec<String>>();

        //Find variants from bubbles
        let vcf_list = detect_all_variants(
            &path_to_steps_map,
            &possible_bubbles_list,
            &graph,
            &node_id_to_path_and_pos_map,
            false,
            100,
            &paths_list,
        );

        (vcf_list, g_bfs, possible_bubbles_list)
    }

    #[test]
    fn test_variant_detection() {
        let graph = read_test_gfa();
        let variants_found = run_whole_script(graph).0;

        //Check that all variants have been found
        assert_eq!(variants_found.len(), 21);

        //Check one variant per type

        //x	9	.	G	A	.	.	TYPE=snv	GT	0|1
        let snv = Variant {
            chromosome: "x".to_string(),
            position: 9,
            id: None,
            reference: "G".to_string(),
            alternate: Some("A".to_string()),
            quality: None,
            filter: None,
            info: Some("TYPE=snv".to_string()),
            format: Some("GT".to_string()),
            sample_name: Some("0|1".to_string()),
        };
        assert!(variants_found.contains(&snv));

        //x	18	.	T	TAA	.	.	TYPE=ins	GT	0|1
        let ins = Variant {
            chromosome: "x".to_string(),
            position: 18,
            id: None,
            reference: "T".to_string(),
            alternate: Some("TAA".to_string()),
            quality: None,
            filter: None,
            info: Some("TYPE=ins".to_string()),
            format: Some("GT".to_string()),
            sample_name: Some("0|1".to_string()),
        };
        assert!(variants_found.contains(&ins));

        //y	18	.	TAA	T	.	.	TYPE=del	GT	0|1
        let del = Variant {
            chromosome: "y".to_string(),
            position: 18,
            id: None,
            reference: "TAA".to_string(),
            alternate: Some("T".to_string()),
            quality: None,
            filter: None,
            info: Some("TYPE=del".to_string()),
            format: Some("GT".to_string()),
            sample_name: Some("0|1".to_string()),
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
        graph.create_edge(&Edge(h1, h2));
        graph.create_edge(&Edge(h2, h3));
        graph.create_edge(&Edge(h3, h4));
        //Loop
        graph.create_edge(&Edge(h3, h5));
        graph.create_edge(&Edge(h5, h3));

        let paths = find_all_paths_between(&graph, &h1.id(), &h4.id(), 100);

        assert!(paths.len() == 1);
        assert!(paths.contains(&vec![h1.id(), h2.id(), h3.id(), h4.id()]));
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
        graph.create_edge(&Edge(h1, h2));
        //Path 2
        graph.create_edge(&Edge(h1, h3));
        graph.create_edge(&Edge(h3, h2));
        //Path 3
        graph.create_edge(&Edge(h3, h4));
        graph.create_edge(&Edge(h4, h2));

        let paths = find_all_paths_between(&graph, &h1.id(), &h2.id(), 100);

        assert!(paths.len() == 3);
        assert!(paths.contains(&vec![h1.id(), h2.id()]));
        assert!(paths.contains(&vec![h1.id(), h3.id(), h2.id()]));
        assert!(paths.contains(&vec![h1.id(), h3.id(), h4.id(), h2.id()]));
    }

    #[test]
    fn find_all_paths_3() {
        let mut graph = HashGraph::new();

        //Add nodes
        let h1 = graph.append_handle("A");
        let h2 = graph.append_handle("T");
        let h3 = graph.append_handle("C");
        let h4 = graph.append_handle("G");
        let h5 = graph.append_handle("AT");

        //Add edges
        //Path 1
        graph.create_edge(&Edge(h1, h2));
        graph.create_edge(&Edge(h2, h4));
        //Path 2
        graph.create_edge(&Edge(h1, h3));
        graph.create_edge(&Edge(h3, h4));
        
        graph.create_edge(&Edge(h4, h5));

        let paths = find_all_paths_between(&graph, &h1.id(), &h5.id(), 100);

        assert!(paths.len() == 2);
        assert!(paths.contains(&vec![h1.id(), h2.id(), h4.id(), h5.id()]));
        assert!(paths.contains(&vec![h1.id(), h3.id(), h4.id(), h5.id()]));
    }

    #[test]
    fn find_all_paths_4() {
        let mut graph = HashGraph::new();

        //Add nodes
        let h1 = graph.append_handle("A");
        let h2 = graph.append_handle("T");
        let h3 = graph.append_handle("C");
        let h4 = graph.append_handle("G");
        let h5 = graph.append_handle("AT");
        let h6 = graph.append_handle("GT");

        //Add edges
        //Path 1
        graph.create_edge(&Edge(h1, h2));
        graph.create_edge(&Edge(h2, h3));
        graph.create_edge(&Edge(h3, h5));
        //Path 2
        graph.create_edge(&Edge(h1, h4));
        graph.create_edge(&Edge(h4, h5));
        
        graph.create_edge(&Edge(h5, h6));

        let paths = find_all_paths_between(&graph, &h1.id(), &h6.id(), 100);

        assert!(paths.len() == 2);
        assert!(paths.contains(&vec![h1.id(), h2.id(), h3.id(), h5.id(), h6.id()]));
        assert!(paths.contains(&vec![h1.id(), h4.id(), h5.id(), h6.id()]));
    }

    #[test]
    fn find_all_paths_4_alt() {
        let mut graph = HashGraph::new();

        //Add nodes
        let h1 = graph.append_handle("A");
        let h2 = graph.append_handle("T");
        let h3 = graph.append_handle("C");
        let h4 = graph.append_handle("G");
        let h5 = graph.append_handle("AT");
        let h6 = graph.append_handle("GT");

        //Add edges
        //Path 1
        graph.create_edge(&Edge(h1, h2));
        graph.create_edge(&Edge(h2, h3));
        graph.create_edge(&Edge(h3, h5));
        //Path 2
        graph.create_edge(&Edge(h1, h4));
        graph.create_edge(&Edge(h4, h5));
        
        graph.create_edge(&Edge(h5, h6));

        let paths = find_all_paths_between(&graph, &h1.id(), &h6.id(), 100);

        assert!(paths.len() == 2);
        assert!(paths.contains(&vec![h1.id(), h2.id(), h3.id(), h5.id(), h6.id()]));
        assert!(paths.contains(&vec![h1.id(), h4.id(), h5.id(), h6.id()]));
    }

    #[test]
    fn find_all_paths_5() {
        let mut graph = HashGraph::new();

        //Add nodes
        let h1 = graph.append_handle("A");
        let h2 = graph.append_handle("T");
        let h3 = graph.append_handle("C");
        let h4 = graph.append_handle("G");
        let h5 = graph.append_handle("AT");
        let h6 = graph.append_handle("GT");

        //Add edges
        //Path 1
        graph.create_edge(&Edge(h1, h2));

        graph.create_edge(&Edge(h2, h3));
        graph.create_edge(&Edge(h2, h4));

        graph.create_edge(&Edge(h3, h5));
        graph.create_edge(&Edge(h4, h5));
        graph.create_edge(&Edge(h5, h6));

        //Path 2
        graph.create_edge(&Edge(h1, h6));

        let paths = find_all_paths_between(&graph, &h1.id(), &h6.id(), 100);

        assert!(paths.len() == 3);
        assert!(paths.contains(&vec![h1.id(), h2.id(), h3.id(), h5.id(), h6.id()]));
        assert!(paths.contains(&vec![h1.id(), h2.id(), h4.id(), h5.id(), h6.id()]));
        assert!(paths.contains(&vec![h1.id(), h6.id()]));  
    }

    #[test]
    fn find_all_paths_5_alt() {
        let mut graph = HashGraph::new();

        //Add nodes
        let h1 = graph.append_handle("A");
        let h2 = graph.append_handle("T");
        let h3 = graph.append_handle("C");
        let h4 = graph.append_handle("G");
        let h5 = graph.append_handle("AT");
        let h6 = graph.append_handle("GT");
        let h7 = graph.append_handle("CG");

        //Add edges
        //Path 1
        graph.create_edge(&Edge(h1, h2));

        graph.create_edge(&Edge(h2, h3));
        graph.create_edge(&Edge(h2, h4));

        graph.create_edge(&Edge(h3, h5));
        graph.create_edge(&Edge(h4, h5));
        graph.create_edge(&Edge(h5, h6));
        graph.create_edge(&Edge(h6, h7));

        //Path 2
        graph.create_edge(&Edge(h1, h7));

        let paths = find_all_paths_between(&graph, &h1.id(), &h7.id(), 100);

        assert!(paths.len() == 3, "There are {} paths insted: {:#?}", paths.len(), paths);
        assert!(paths.contains(&vec![h1.id(), h2.id(), h3.id(), h5.id(), h6.id(), h7.id()]));
        assert!(paths.contains(&vec![h1.id(), h2.id(), h4.id(), h5.id(), h6.id(), h7.id()]));
        assert!(paths.contains(&vec![h1.id(), h7.id()]));  
    }

    #[test]
    fn whole_script_1() {
        // Test the script in a simple linear path
        
        let mut graph = HashGraph::new();

        //Add nodes
        let h1 = graph.append_handle("A");
        let h2 = graph.append_handle("T");
        let h3 = graph.append_handle("C");
        let h4 = graph.append_handle("G");

        //Add edges
        graph.create_edge(&Edge(h1, h2));
        graph.create_edge(&Edge(h2, h3));
        graph.create_edge(&Edge(h3, h4));

        let (variants, g_bfs, bubbles) = run_whole_script(graph);

        // for h in handles_iter(&g_bfs) {
        //     display_node_edges(&g_bfs, &h);
        // }

        assert!(g_bfs.has_edge(h1, h2));
        assert!(g_bfs.has_edge(h2, h3));
        assert!(g_bfs.has_edge(h3, h4));

        assert!(variants.len() == 0);
        assert!(bubbles.len() == 0);
    }

    #[test]
    fn whole_script_2() {
        // Test the script with a simple bubble
        
        let mut graph = HashGraph::new();

        // Add nodes
        let h1 = graph.append_handle("A");
        let h2 = graph.append_handle("T");
        let h3 = graph.append_handle("C");
        let h4 = graph.append_handle("G");

        // Add edges
        graph.create_edge(&Edge(h1, h2));
        graph.create_edge(&Edge(h1, h3));
        graph.create_edge(&Edge(h2, h4));
        graph.create_edge(&Edge(h3, h4));

        // Add a reference path
        let path1 = graph.create_path_handle("x", false);
        graph.append_step(&path1, h1);
        graph.append_step(&path1, h2);
        graph.append_step(&path1, h4);
        // println!("Paths are: {:#?}",graph.paths);

        let (variants, g_bfs, bubbles) = run_whole_script(graph);

        // g_bfs
        assert!(g_bfs.has_edge(h1, h2));
        assert!(g_bfs.has_edge(h1, h3));
        assert!(g_bfs.has_edge(h3, h4) || g_bfs.has_edge(h2, h4));
        assert!(!(g_bfs.has_edge(h3, h4) && g_bfs.has_edge(h2, h4)));

        // bubbles
        assert!(bubbles.len() == 1);
        let bubble_1 = Bubble {
            start : h1.id(),
            end : h4.id()
        };
        assert!(bubbles.contains(&bubble_1));
        
        // variants
        assert!(variants.len() == 1);
        let variant_1 = Variant {
            chromosome: "x".to_string(),
            position: 2,
            id: None,
            reference: "T".to_string(),
            alternate: Some(
                "C".to_string(),
            ),
            quality: None,
            filter: None,
            info: Some(
                "TYPE=snv".to_string(),
            ),
            format: Some(
                "GT".to_string(),
            ),
            sample_name: Some(
                "0|1".to_string(),
            ),
        };
        assert!(variants.contains(&variant_1));
    }

    #[test]
    fn whole_script_3() {
        // Test the script with a nested bubble
        
        let mut graph = HashGraph::new();

        // Add nodes
        let h1 = graph.append_handle("A");
        let h2 = graph.append_handle("T");
        let h3 = graph.append_handle("C");
        let h4 = graph.append_handle("G");
        let h5 = graph.append_handle("AT");
        let h6 = graph.append_handle("GT");

        // Add edges
        graph.create_edge(&Edge(h1, h2));
        graph.create_edge(&Edge(h1, h6));

        graph.create_edge(&Edge(h2, h3));
        graph.create_edge(&Edge(h2, h4));

        graph.create_edge(&Edge(h3, h5));
        graph.create_edge(&Edge(h4, h5));
        graph.create_edge(&Edge(h5, h6));

        // Add a reference path
        let path1 = graph.create_path_handle("y", false);
        graph.append_step(&path1, h1);
        graph.append_step(&path1, h6);

        let (variants, g_bfs, bubbles) = run_whole_script(graph);

        // g_bfs
        // External bubble
        assert!(g_bfs.has_edge(h1, h2));
        assert!(g_bfs.has_edge(h1, h6));
        
        // Internal bubble
        assert!(g_bfs.has_edge(h2, h3));
        assert!(g_bfs.has_edge(h2, h4));

        // Close internal bubble
        assert!(g_bfs.has_edge(h3, h5) || g_bfs.has_edge(h4, h5));
        assert!(!(g_bfs.has_edge(h3, h5) && g_bfs.has_edge(h4, h5)));

        //Close external bubble
        assert!(!(g_bfs.has_edge(h5, h6)));
        

        // bubbles
        assert!(bubbles.len() == 2, "Bubbles has len {} and is {:#?}", bubbles.len(), bubbles);
        let bubble_1 = Bubble {
            start : h1.id(),
            end : h6.id()
        };
        let bubble_2 = Bubble {
            start : h2.id(),
            end : h5.id()
        };
        assert!(bubbles.contains(&bubble_1));
        assert!(bubbles.contains(&bubble_2));
        
        // variants
        assert!(variants.len() == 1);
        let variant_1 = Variant {
            chromosome: "x".to_string(),
            position: 2,
            id: None,
            reference: "T".to_string(),
            alternate: Some(
                "C".to_string(),
            ),
            quality: None,
            filter: None,
            info: Some(
                "TYPE=snv".to_string(),
            ),
            format: Some(
                "GT".to_string(),
            ),
            sample_name: Some(
                "0|1".to_string(),
            ),
        };
        assert!(variants.contains(&variant_1));
    }


}
