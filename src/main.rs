use std::path::PathBuf;
fn main() {
    println!("Hello, world!");
    let gfa = gfa::parser::parse_gfa(&PathBuf::from("./input/lil.gfa"));
    
    match gfa {
        None => panic!("Error parsing GFA file"),
        Some(g) => {
            let num_segs = g.segments.len();
            let num_links = g.links.len();
            let num_paths = g.paths.len();
            let num_conts = g.containments.len();

            println!("num_segs: {}",num_segs);
            println!("num_links: {}",num_links);
            println!("num_paths: {}",num_paths);
            println!("num_conts: {}",num_conts);
        }
    }
}
