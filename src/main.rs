use std::path::PathBuf;
use petgraph::Graph;
use petgraph::dot::Dot;

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

            let mut deps = Graph::<&str,()>::new();
            
            println!("Segments");
            for seg in &g.segments {
                println!("{},{}",seg.name,seg.sequence);
            }

            println!("Links");
            for link in &g.links {
                println!("{},{},{},{},{}",link.from_segment, link.from_orient.as_bool(), link.to_segment, link.to_orient.as_bool(), link.overlap);
            }

            println!("Paths");
            for path in &g.paths {
                println!("{},{:?},{:?}", path.path_name, path.segment_names, path.overlaps);
            }
            
            // println("Containments");
            // for c in &g.containments {
            //     println!("")
            // }


            //println!("{}", Dot::new(&deps));

        }
    }
}
