use std::path::PathBuf;
use handlegraph::graph::HashGraph;
use std::fs::File;
use std::io::Write;
extern crate clap;
use clap::{Arg, App};

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

    let matches = App::new("rs-GFAtoVCF")
                          .version("1.0")
                          .author("Francesco Porto <francesco.porto97@gmail.com>")
                          .about("Converts GFA to VCF")
                          .arg(Arg::with_name("input")
                               .short("i")
                               .long("input")
                               .value_name("FILE")
                               .help("Sets the input file to use")
                               .required(true)
                               .takes_value(true))
                          .arg(Arg::with_name("output")
                               .short("o")
                               .long("output")
                               .value_name("FILE")
                               .help("Sets the output file to use")
                               .required(true)
                               .takes_value(true))
                          .get_matches();

    let in_path = matches.value_of("input").expect("Missing parameter --input");
    let out_path = matches.value_of("output").expect("Missing parameter --output");

    //let in_path = "./input/samplePath3.gfa";
    //let out_path = "./input/samplePath3.vcf";

    //let gfa = gfa::parser::parse_gfa(&PathBuf::from(in_path)).unwrap();

    if let Some(gfa) = gfa::parser::parse_gfa(&PathBuf::from(in_path)) {
        let graph = HashGraph::from_gfa(&gfa);
        //println!("{:?}",graph);
        let variations = find_variants(&graph);
        write_to_file(&PathBuf::from(out_path), &variations).unwrap();
    } else {
        panic!("Couldn't parse test GFA file!");
    }

}

fn find_variants(graph: &HashGraph) -> Vec<Variant> {
    // TODO: function
    let mut variations: Vec<Variant> = Vec::new();

    //let node_ids: Vec<_> = graph.graph.keys().collect();
            //println!("Node IDs:");

    variations
}

fn write_to_file(path: &PathBuf, variations: &Vec<Variant>) -> std::io::Result<()> {
    let mut file = File::create(path).expect(&format!("Error creating file {:?}", path));

    //TODO: header, ...

    let fields = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    file.write(fields.as_bytes()).expect("Error writing header");
    for var in variations {
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
