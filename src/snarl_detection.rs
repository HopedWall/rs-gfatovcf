use handlegraph::hashgraph::HashGraph;
use handlegraph::handle::{Handle,Direction};
use handlegraph::handlegraph::{handles_iter,handle_edges_iter};
use std::collections::HashMap;
use std::collections::HashSet;
use std::collections::VecDeque;

type HandleSet = HashSet<Handle>;
type Handle2Component = HashMap<Handle, usize>;

/// Step 1) Get undirected adjacency connected components of VG *sides*
fn compute_side_components(graph : &HashGraph, components : &mut Vec<HandleSet>, handle_to_component : &mut Handle2Component) {

    let mut add_handle = |handle| {
         if !handle_to_component.contains_key(handle) {
            
            // Make a new component
            let component_id = components.len();
            //emplace back?
            let component = components.last_mut().unwrap();

            // Mark the current node as being on this component
            *handle_to_component.get_mut(handle).unwrap() = component_id;
            component.insert(*handle);

            //Do a BFS
            let mut q: VecDeque<Handle> = VecDeque::new();
            q.push_back(*handle);

            while let Some(next_handle) = q.pop_front() {
                for traversed_handle in handle_edges_iter(graph, next_handle, Direction::Right) {
                    let adjacent_side = traversed_handle.flip();

                    if !handle_to_component.contains_key(&adjacent_side) {
                        *handle_to_component.get_mut(&adjacent_side).unwrap() = component_id;
                        component.insert(adjacent_side);
                        q.push_back(adjacent_side);
                    }
                }
            }
         }
    };

    
    // Collect handles early otherwise they do not live long enough
    let graph_handles : Vec<Handle> = handles_iter(graph).collect();
    let flipped_handles : Vec<Handle> = graph_handles.iter().map(|x| x.flip()).collect();

    for i in 0..graph_handles.len() {
        add_handle(&graph_handles[i]);
        add_handle(&flipped_handles[i]);
    }

}



