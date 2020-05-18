use std::path::PathBuf;
use handlegraph::graph::HashGraph;
use std::fs::File;
use std::io::Write;

struct Variant {
    chromosome: String,
    position: i32,
    id: String,
    reference: String,
    alternate: String,
    quality: i32,
    filter: String,
    info: String,
}

fn main() {
    
    println!("Insert path to file");
    
    // let mut path = String::new();
    // std::io::stdin().read_line(&mut path).expect("Error");
    // path.pop(); //Remove trailing /n

    let path = "./input/samplePath3.gfa";
    let out = "./output/samplePath3.vcf";

    let gfa = gfa::parser::parse_gfa(&PathBuf::from(path)).unwrap();

    let graph = HashGraph::from_gfa(&gfa);
    
    println!("{:?}",graph);

    // let mut edges: Vec<(&str, &str)> = Vec::new();
    // for link in &gfa.links {
    //     edges.push((&link.from_segment, &link.to_segment))
    // }
    // let g = DiGraphMap::<&str, &str>::from_edges(edges);
    // let x = Dot::with_config(&g, &[Config::EdgeIndexLabel]);
    // println!("{}", x);

    let variations = find_variants(graph);
    write_to_file(&PathBuf::from(out), variations).unwrap();

}

fn find_variants(graph: HashGraph) -> Vec<Variant> {
    // TODO: function
    let mut variations: Vec<Variant> = Vec::new();
    variations
}

fn write_to_file(path: &PathBuf, variations: Vec<Variant>) -> std::io::Result<()> {
    let mut file = File::create(path).expect(&format!("Error creating file {:?}", path));

    //TODO: header, ...

    let fields = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    file.write(fields.as_bytes()).expect("Error writing header");
    for var in &variations {
        let to_write = format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                                var.chromosome, 
                                var.position, 
                                var.id, 
                                var.reference, 
                                var.alternate, 
                                var.quality, 
                                var.filter, 
                                var.info);
        file.write(to_write.as_bytes()).expect("Error writing variant");
    }
    Ok(())
}
