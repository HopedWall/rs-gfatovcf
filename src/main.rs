use std::path::PathBuf;
use petgraph::dot::Dot;
// for the graph
use petgraph::dot::Config;
use petgraph::graphmap::DiGraphMap;

fn main() {
    
    println!("Insert path to file");
    
    // let mut path = String::new();
    // std::io::stdin().read_line(&mut path).expect("Error");
    // path.pop(); //Remove trailing /n

    let path = "./input/samplePath3.gfa";

    let gfa = gfa::parser::parse_gfa(&PathBuf::from(path)).unwrap();

    let mut edges: Vec<(&str, &str)> = Vec::new();
    for link in &gfa.links {
        edges.push((&link.from_segment, &link.to_segment))
    }
    let g = DiGraphMap::<&str, &str>::from_edges(edges);
    let x = Dot::with_config(&g, &[Config::EdgeIndexLabel]);
    println!("{}", x);

}
