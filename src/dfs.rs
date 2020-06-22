/// A few ideas for dfs-based detection
/// CURRENTLY NOT IN USE!!!

/// OLD PYTHON DFS
/// Computes the DFS of a given HashGraph
// fn dfs(g : &HashGraph, g_dfs : &mut HashGraph, node_id : &NodeId) {
//     let current_node = g.get_handle(*node_id, false);
//     //let sequence_node = g.get_sequence(&current_node);

//     g.follow_edges(
//         &current_node, 
//         Direction::Right,  //What should go here?
//         |neighbor| {
//             create_edge_and_so_on(
//                 &g, g_dfs, 
//                 &current_node, neighbor,
//                 &g.get_id(neighbor));
//                 true
//         });
// }
/// Adds edges to DFS
// fn create_edge_and_so_on(g : &HashGraph, g_dfs : &mut HashGraph, handle1 : &Handle, handle2 : &Handle, neighbor_node_id : &NodeId) {
//     let handle1_id = g.get_id(handle1);
//     let handle2_id = g.get_id(handle2); 
//     //println!("Create edge from {} to {}",handle1_id,handle2_id);

//     if !g_dfs.has_node(handle2_id) {
//         //println!("g_dfs does not have node: {}",handle2_id);
//         dfs(g,g_dfs,neighbor_node_id);
        
//         if !g_dfs.has_node(handle1_id) {
//             //println!("Add node {}",handle1_id);
//             g_dfs.create_handle(
//                 g.get_sequence(handle1), 
//                 handle1_id);
//         }
    
//         if !g_dfs.has_node(handle2_id) {
//             //println!("Add node {}",handle2_id);
//             g_dfs.create_handle(
//                 g.get_sequence(handle2), 
//                 handle2_id);
//         }
    
//         g_dfs.create_edge(handle1, handle2)
//     }
// }

//Computes a different dfs
// fn new_dfs(g : &HashGraph, g_dfs : &HashGraph, node_id : &NodeId) {
//     let current_handle = g.get_handle(*node_id, false);

//     g.follow_edges(
//         &current_handle,
//         Direction::Right,
//         |neighbor| {
//             if !g_dfs.has_node(g.get_id(neighbor)) {
//                 new_dfs_support(g, g_dfs, &neighbor);
//             }         
//             true
//         });
// }
// fn new_dfs_support(g : &HashGraph, g_dfs : &HashGraph, node_handle : &Handle) {
//}

// Find all paths not in g_dfs but in dfs
// fn find_bubbles_dfs(g : &HashGraph, g_dfs : &HashGraph, node_id : &NodeId, bubbles_vector : &mut Vec<(NodeId, NodeId)>) {
//     let current_handle = g.get_handle(*node_id, false);

//     g.follow_edges(
//         &current_handle, 
//         Direction::Right,
//         |neighbor| {
//             if !g_dfs.has_edge(&current_handle, neighbor) {
//                 let possible_bubble : (NodeId, NodeId) = (*node_id,NodeId::from(0)); 
//                 find_bubbles_dfs_support(&g, &g_dfs, &neighbor, bubbles_vector, possible_bubble);
//             }
//             //find_bubbles(g, g_dfs, &g.get_id(neighbor), bubbles_vector);
//             true
//     });    
// }
// // Once a bubble opening has been found, continue exploring the graph
// fn find_bubbles_dfs_support(g : &HashGraph, g_dfs : &HashGraph, node_handle : &Handle, 
//                         bubbles_vector : &mut Vec<(NodeId, NodeId)>,
//                         mut current_bubble : (NodeId, NodeId)) {
        
//     g.follow_edges(
//         node_handle, 
//         Direction::Right, 
//         |neighbor| {
//             if g_dfs.has_node(g.get_id(neighbor)) {
//                 // Set bubble end
//                 current_bubble.1 = g.get_id(neighbor);
//                 // Add bubble to vector
//                 if !bubbles_vector.contains(&current_bubble) {
//                     bubbles_vector.push(current_bubble);
//                 }
//                 //Continue searching new bubble
//                 find_bubbles_dfs(g, g_dfs, &g.get_id(neighbor), bubbles_vector);
//             } else {
//                 find_bubbles_dfs_support(g, g_dfs, neighbor, 
//                     bubbles_vector, current_bubble);
//             }
//             true
//         });
// }

/// TESTS used in dfs-version
// mod tests {
//     use handlegraph::graph::HashGraph;
//     use handlegraph::handlegraph::HandleGraph;
//     use super::*;

//     //Used in other tests
//     fn read_test_gfa() -> HashGraph {
//         use gfa::parser::parse_gfa;

//         HashGraph::from_gfa(&parse_gfa(&PathBuf::from("./input/samplePath3.gfa")).unwrap())
//     }

//     #[test]
//     fn test_dfs_1() {
//         let mut graph = HashGraph::new();
//         let h1 = graph.append_handle("A");
//         let h2 = graph.append_handle("CG");

//         graph.create_edge(&h1, &h2);

//         let mut g_dfs = HashGraph::new();
//         dfs(&graph,&mut g_dfs,&h1.id());

//         // g_dfs.for_each_handle(|h| {
//         //     display_node_edges(&g_dfs, &h);
//         //     true
//         // });
        
//     }

//     #[test]
//     fn test_dfs_2() {

//         // Test a pattern that provides a different dfs_tree than Flavia's

//         /*
//         edges
//         1  ----> 3 
//           \-> 2 /
//         */

//         let mut graph = HashGraph::new();
//         let h1 = graph.append_handle("A");
//         let h2 = graph.append_handle("CG");
//         let h3 = graph.append_handle("T");


//         graph.create_edge(&h1, &h2);
//         graph.create_edge(&h1, &h3);
//         graph.create_edge(&h2, &h3);
        

//         let mut g_dfs = HashGraph::new();
//         dfs(&graph,&mut g_dfs,&h1.id());

//         // g_dfs.for_each_handle(|h| {
//         //     display_node_edges(&g_dfs, &h);
//         //     true
//         // });

//         //println!("{:?}",graph);
//         //println!("{:?}",g_dfs);
//     }

//     #[test]
//     fn test_dfs_3() {

//         // Test a pattern that provides the same dfs_tree as Flavia's

//         /*
//         edges
        
//           /-> 3 -\ 
//         1         4
//           \-> 2 -/
//         */

//         let mut graph = HashGraph::new();
//         let h1 = graph.append_handle("A");
//         let h2 = graph.append_handle("CG");
//         let h3 = graph.append_handle("T");
//         let h4 = graph.append_handle("AC");

//         graph.create_edge(&h1, &h2);
//         graph.create_edge(&h1, &h3);
//         graph.create_edge(&h2, &h4);
//         graph.create_edge(&h3, &h4);

//         let mut g_dfs = HashGraph::new();
//         dfs(&graph,&mut g_dfs,&h1.id());

//         //println!("{:?}",graph);

//         // g_dfs.for_each_handle(|h| {
//         //     display_node_edges(&g_dfs, &h);
//         //     true
//         // });

//         //println!("{:?}",g_dfs);
//     }
// }