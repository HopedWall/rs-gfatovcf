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
use log::info;
use std::collections::HashSet;
use std::collections::VecDeque;

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

/// Performs a BFS over a given HashGraph, and returns the resulting bfs-tree
/// See section 'Example' from https://en.wikipedia.org/wiki/Breadth-first_search
/// for more details
pub fn bfs(g: &HashGraph, node_id: &NodeId) -> HashGraph {
    // The resulting tree will still be stored inside a HashGraph for the
    // sake of simplicity; note that this is in fact a tree, i.e.
    // it has a root, there are no loops, etc.
    let mut g_bfs = HashGraph::new();

    // Create queue
    // NOTE: this is a Queue based implementation, this was done
    // in order not to get a stack overflow (the previous recursion-based
    // version was often experiencing this kind of issue)
    let mut q: VecDeque<NodeId> = VecDeque::new();

    // Insert first value
    q.push_back(*node_id);

    while !q.is_empty() {
        info!("Queue is {:#?}", q);

        let curr_node = q.pop_front().unwrap();
        let current_handle = Handle::pack(curr_node, false);

        info!("Curr node is {:#?}", curr_node);

        // Check if curr_node is already in g_bfs
        if !g_bfs.has_node(curr_node) {
            g_bfs.create_handle(g.sequence(current_handle), curr_node);
        }

        for neighbor in handle_edges_iter(g, current_handle, Direction::Right) {
            if !g_bfs.has_node(neighbor.id()) {
                // Create handle in g_bfs
                g_bfs.create_handle(g.sequence(neighbor), neighbor.id());

                // Add neighbor id to queue
                q.push_back(neighbor.id());

                // Create edge from curr_handle to new node in g_bfs
                let edge = Edge::edge_handle(current_handle, neighbor);
                g_bfs.create_edge(&edge);

                // Add new node to queue
                q.push_back(neighbor.id());
            }
        }
    }

    g_bfs
}

/// Prints a node of a given HashGraph, and all of its edges
pub fn display_node_edges(g_dfs: &HashGraph, h: &Handle) {
    println!("node {}", h.id());

    // Prints all the outgoing edges from the current node
    for n in handle_edges_iter(g_dfs, *h, Direction::Right) {
        println!("{} --> {}", h.id(), n.id());
    }
}

/// Finds the distance of each node from a given root
pub fn bfs_distances(
    g: &HashGraph,
    starting_node_id: &NodeId,
) -> (BTreeMap<NodeId, u64>, Vec<NodeId>) {
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
        info!("Queue is {:#?}", q);
        info!("Curr node id is {:#?}", current_node_id);

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

/// Returns a map where, for each distance from the root (i.e. an integer between 0 and the maximum distance),
/// contains the number of nodes at that distance
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

/// Detects bubbles in a variation graph. A bubble consists of multiple directed unipaths
/// (a path in which all internal vertices are of degree 2) between two vertices. In the
/// context of this program, bubbles are represented as pairs (Start_Node_Id, End_Node_Id)
pub fn detect_bubbles(
    distances_map: &BTreeMap<NodeId, u64>,
    ordered_node_id_list: &[NodeId],
    dist_to_num_nodes: &BTreeMap<u64, usize>,
) -> Vec<(NodeId, NodeId)> {
    let mut possible_bubbles_list: Vec<(NodeId, NodeId)> = Vec::new();
    let mut open_bubble = false;

    let mut curr_bubble: (NodeId, NodeId) = (NodeId::from(0), NodeId::from(0));
    for node_id in ordered_node_id_list {
        // Get distance of current NodeId from root in g_bfs
        let node_distance = distances_map[&node_id];

        // If there is only 1 node at that distance
        if dist_to_num_nodes[&node_distance] == 1 {
            // And there are multiple nodes at distance+1 -> open bubble
            // Note: node_distance could go out of bounds
            if ((node_distance + 1) as usize) < dist_to_num_nodes.len()
                && dist_to_num_nodes[&(node_distance + 1)] > 1
            {
                // Close current bubble if one is already open
                if open_bubble {
                    curr_bubble.1 = *node_id;
                    possible_bubbles_list.push(curr_bubble);
                    curr_bubble = (NodeId::from(0), NodeId::from(0));
                }

                // Start new bubble
                curr_bubble.0 = *node_id;
                open_bubble = true;
            } else {
                // If a bubble is open
                if open_bubble {
                    // Close bubble
                    curr_bubble.1 = *node_id;
                    possible_bubbles_list.push(curr_bubble);

                    // Reset curr_bubble for future bubbles
                    curr_bubble = (NodeId::from(0), NodeId::from(0));
                    open_bubble = false;
                }
            }
        }
    }

    possible_bubbles_list
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
mod tests {
    use super::*;
    use handlegraph::handlegraph::HandleGraph;
    use handlegraph::hashgraph::HashGraph;
    use std::path::PathBuf;


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
        graph.create_edge(&Edge(h1, h2));
        graph.create_edge(&Edge(h2, h3));
        graph.create_edge(&Edge(h3, h4));
        //Loop
        graph.create_edge(&Edge(h3, h5));
        graph.create_edge(&Edge(h5, h3));
    
        let g_bfs = bfs(&graph, &NodeId::from(1));
    
        //Check g_bfs does not contain the loop
        assert!(g_bfs.has_edge(h3, h5));
        assert!(!g_bfs.has_edge(h5, h3));
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
        graph.create_edge(&Edge(h1, h2));
        //Path 2
        graph.create_edge(&Edge(h1, h3));
        graph.create_edge(&Edge(h3, h2));
        //Path 3
        graph.create_edge(&Edge(h3, h4));
        graph.create_edge(&Edge(h4, h2));
    
        let g_bfs = bfs(&graph, &NodeId::from(1));
    
        assert!(g_bfs.has_edge(h1, h2));
        assert!(g_bfs.has_edge(h1, h3));
        assert!(g_bfs.has_edge(h3, h4));
    
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
        graph.create_edge(&Edge(h1, h2));
        graph.create_edge(&Edge(h2, h3));
        graph.create_edge(&Edge(h3, h4));
        //Loop
        graph.create_edge(&Edge(h3, h5));
        graph.create_edge(&Edge(h5, h3));
    
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
        graph.create_edge(&Edge(h1, h2));
        //Path 2
        graph.create_edge(&Edge(h1, h3));
        graph.create_edge(&Edge(h3, h2));
        //Path 3
        graph.create_edge(&Edge(h3, h4));
        graph.create_edge(&Edge(h4, h2));
    
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
    
        assert!(possible_bubbles_list.contains(&(NodeId::from(1), NodeId::from(4))));
    }

}
