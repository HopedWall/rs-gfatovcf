use std::path::PathBuf;
use petgraph::Graph;
use petgraph::dot::Dot;
// for the graph
use petgraph::dot::Config;
use petgraph::graphmap::DiGraphMap;

fn main() {
    
    println!("Insert path to file");
    
    // let mut path = String::new();
    // std::io::stdin().read_line(&mut path).expect("Error");
    // path.pop(); //Remove trailing /n

    let mut path = "./input/samplePath3.gfa";

    let gfa = gfa::parser::parse_gfa(&PathBuf::from(path));
    
    match gfa {
        None => panic!("Error parsing GFA file"),
        Some(g) => {

            //let mut graph = Graph::<(&str, &str),_>::new();
            //let mut graph = Graph::<&str,&str>::new();

            println!("Segments");
            for seg in &g.segments {
                //println!("{},{}",seg.name,seg.sequence);
                //let t = (&seg.name, &seg.sequence);
                //graph.add_node(&seg.name);
                
            }

            println!("Links");
            let mut edges: Vec<(&str, &str)> = Vec::new();
            for link in &g.links {
                //println!("{},{},{},{},{}",link.from_segment, link.from_orient.as_bool(), link.to_segment, link.to_orient.as_bool(), link.overlap);
                //graph.add_edge(&link.from_segment, &link.to_segment, "0");
                edges.push((&link.from_segment, &link.to_segment))
            }
            let g = DiGraphMap::<&str, &str>::from_edges(edges);
            let x = Dot::with_config(&g, &[Config::EdgeIndexLabel]);
            println!("{}", x);

            // println!("Paths");
            // for path in &g.paths {
            //     println!("{},{:?},{:?}", path.path_name, path.segment_names, path.overlaps);
            // }

            //println!("{}", Dot::new(&graph));

        }
    }
}
