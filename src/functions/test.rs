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

fn run_whole_script(graph: HashGraph) -> Vec<Variant> {
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
        100
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
