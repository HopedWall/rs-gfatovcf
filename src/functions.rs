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
extern crate chrono;
use std::cmp;
use std::cmp::Ordering;
use std::collections::HashSet;
use std::collections::VecDeque;
//use std::io;

/// A struct that holds Variants, as defined in the VCF format
#[derive(PartialEq)]
pub struct Variant {
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

impl Variant {
    pub fn to_vcf_string(&self) -> String {
        let vcf_string = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            self.chromosome,
            self.position,
            self.id,
            self.reference,
            self.alternate,
            self.quality,
            self.filter,
            self.info,
            self.format,
            self.sample_name
        );
        
        vcf_string
    }

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
pub fn paths_to_steps(graph: &HashGraph) -> HashMap<String, Vec<String>> {
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
pub fn bfs(g: &HashGraph, node_id: &NodeId) -> HashGraph {
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
pub fn display_node_edges(g_dfs: &HashGraph, h: &Handle) {
    println!("node {}", h.id());

    for n in handle_edges_iter(g_dfs, *h, Direction::Right) {
        show_edge(h, &n);
    }
}

/// Finds the distance of each node from a given root
pub fn bfs_distances(g: &HashGraph, starting_node_id: &NodeId) -> (BTreeMap<NodeId, u64>, Vec<NodeId>) {
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
pub fn get_dist_to_num_nodes(distances_map: &BTreeMap<NodeId, u64>) -> BTreeMap<u64, usize> {
    let mut dist_to_num_nodes: BTreeMap<u64, usize> = BTreeMap::new();

    distances_map.values().for_each(|dist| {
        *dist_to_num_nodes.entry(*dist).or_default() += 1;
    });

    dist_to_num_nodes
}

/// Returns paths as sequences
pub fn get_path_to_sequence(
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
pub fn detect_bubbles(
    distances_map: &BTreeMap<NodeId, u64>,
    ordered_node_id_list: &[NodeId],
    dist_to_num_nodes: &BTreeMap<u64, usize>,
) -> Vec<(NodeId, NodeId)> {
    let mut possible_bubbles_list: Vec<(NodeId, NodeId)> = Vec::new();
    let mut open_bubble = false;

    let mut curr_bubble: (NodeId, NodeId) = (NodeId::from(0), NodeId::from(0));
    for node_id in ordered_node_id_list {
        //Get distance of current NodeId from root in g_bfs
        let node_distance = distances_map[&node_id];

        //If there is only 1 node at that distance
        if dist_to_num_nodes[&node_distance] == 1 {
            //And there are multiple nodes at distance+1 -> open bubble
            //Note: node_distance could go out of bounds
            if ((node_distance + 1) as usize) < dist_to_num_nodes.len()
                && dist_to_num_nodes[&(node_distance + 1)] > 1
            {
                //Close current bubble if one is already open
                if open_bubble {
                    curr_bubble.1 = *node_id;
                    possible_bubbles_list.push(curr_bubble);
                    curr_bubble = (NodeId::from(0), NodeId::from(0));
                }

                // Start new bubble
                curr_bubble.0 = *node_id;
                open_bubble = true;
            } else {
                //If a bubble is open
                if open_bubble {
                    //Close bubble
                    curr_bubble.1 = *node_id;
                    possible_bubbles_list.push(curr_bubble);

                    //Reset curr_bubble for future bubbles
                    curr_bubble = (NodeId::from(0), NodeId::from(0));
                    open_bubble = false;
                }
            }
        }
    }

    possible_bubbles_list
}

/// Prints all paths in a given HashGrap, starting from a specific node and ending in another node
pub fn find_all_paths_between(
    g: &HashGraph,
    start_node_id: &NodeId,
    end_node_id: &NodeId,
    max_edges: i32
) -> Vec<Vec<NodeId>> {
    let mut all_paths_list: Vec<Vec<NodeId>> = Vec::new();

    //Put a limit on the maximum amount of edges that can be traversed
    //this should prevent eccessive memory usage
    //let max_edges = 100;
    let mut curr_edges = 0;
    let mut edges_limit_reached = false;

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

            //Break if too many edges have been visited
            curr_edges = curr_edges + 1;
            if curr_edges > max_edges {
                edges_limit_reached = true;
                break;
            }
        }

        if edges_limit_reached {
            break;
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
pub fn detect_all_variants(
    path_to_steps_map: &HashMap<String, Vec<String>>,
    possible_bubbles_list: &[(NodeId, NodeId)],
    graph: &HashGraph,
    node_id_to_path_and_pos_map: &BTreeMap<NodeId, HashMap<String, usize>>,
    verbose: bool,
    max_edges: i32
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
            max_edges
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
    max_edges: i32
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

        //println!("BEFORE FIND ALL PATHS BETWEEN");

        let all_path_list: Vec<Vec<NodeId>> = find_all_paths_between(&graph, start, end, max_edges);

        //println!("AFTER FIND ALL PATHS BETWEEN");

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
pub fn get_node_positions_in_paths(
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

#[cfg(test)]
mod test;