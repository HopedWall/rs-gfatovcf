use std::path::PathBuf;
use handlegraph::graph::HashGraph;
use handlegraph::handlegraph::HandleGraph;
use handlegraph::handle::NodeId;
use gfa::gfa::GFA;
use std::fs::File;
use std::io::Write;
extern crate clap;
use clap::{Arg, App};
extern crate chrono;
use chrono::Utc;

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

    let in_path = matches.value_of("input").expect("Could not parse argument --input");
    let out_path = matches.value_of("output").expect("Could not parse argument --output");

    //let in_path = "./input/samplePath3.gfa";
    //let out_path = "./input/samplePath3.vcf";

    if let Some(gfa) = gfa::parser::parse_gfa(&PathBuf::from(in_path)) {
        let graph = HashGraph::from_gfa(&gfa);

        println!("{:?}",graph);
        let variations = find_variants(&graph);
        write_to_file(&PathBuf::from(out_path), &variations).unwrap();
    } else {
        panic!("Couldn't parse test GFA file!");
    }

}

fn find_variants(graph: &HashGraph) -> Vec<Variant> {
    // TODO: function
    let mut variations: Vec<Variant> = Vec::new();
    
    //is this correct?
    let first_node = graph.get_handle(NodeId::from(1), false);

    //traverse the graph
    //graph.for_each_handle(|h| {
    //   println!();
    //);

    variations
}

fn write_to_file(path: &PathBuf, variations: &Vec<Variant>) -> std::io::Result<()> {
    let mut file = File::create(path).expect(&format!("Error creating file {:?}", path));

    let header = [
        "##fileformat=VCFv4.2",
        &format!("##fileDate={}",Utc::now().format("%Y-%m-%d %H:%M:%S").to_string()),
        "##reference=x.fa",
        "##reference=y.fa",
        "##reference=z.fa",
        "##INFO=<ID=TYPE,Number=A,Type=String,Description=\"Type of each allele (snv, ins, del, mnp, complex)\">",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        &["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"].join("\t"),
    ].join("\n");
    file.write(header.as_bytes()).expect("Error writing header");

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
