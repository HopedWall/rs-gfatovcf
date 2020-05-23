use std::path::PathBuf;
//use handlegraph::graph::HashGraph;
//use handlegraph::handlegraph::HandleGraph;
//use handlegraph::handle::NodeId;
use gfa::gfa::GFA;
use std::fs::File;
use std::io::Write;
extern crate clap;
//use clap::{Arg, App};
extern crate chrono;
use chrono::Utc;
use std::io::{BufReader,BufRead};
use petgraph::graph::Graph;
use petgraph::dot::{Dot, Config};
use petgraph::graphmap::DiGraphMap;
use std::collections::HashMap;

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

    // let matches = App::new("rs-GFAtoVCF")
    //                       .version("1.0")
    //                       .author("Francesco Porto <francesco.porto97@gmail.com>")
    //                       .about("Converts GFA to VCF")
    //                       .arg(Arg::with_name("reference")
    //                            .short("r")
    //                            .long("reference")
    //                            .value_name("FILE")
    //                            .help("Sets the reference file to use")
    //                            .required(true)
    //                            .takes_value(true))
    //                       .arg(Arg::with_name("input")
    //                            .short("i")
    //                            .long("input")
    //                            .value_name("FILE")
    //                            .help("Sets the input file to use")
    //                            .required(true)
    //                            .takes_value(true))
    //                       .arg(Arg::with_name("output")
    //                            .short("o")
    //                            .long("output")
    //                            .value_name("FILE")
    //                            .help("Sets the output file to use")
    //                            .required(true)
    //                            .takes_value(true))
    //                       .get_matches();

    // let in_path = matches.value_of("input").expect("Could not parse argument --input");
    // let out_path = matches.value_of("output").expect("Could not parse argument --output");

    let ref_path = "./input/samplePath3.fa";
    let in_path = "./input/samplePath3.gfa";
    let out_path = "./input/samplePath3.vcf";

    if let Some(reference) = parse_reference(&PathBuf::from(in_path)) {
        
        if let Some(gfa) = gfa::parser::parse_gfa(&PathBuf::from(in_path)) {
            //let graph = HashGraph::from_gfa(&gfa);
            //println!("{:?}",graph);
            
            //Build graph from scratch
            let graph = build_graph(&gfa);

            //Print graph
            let x = Dot::with_config(&graph, &[Config::EdgeIndexLabel]);
            println!("{}", x);

            //Find variations
            let variations = find_variants(&graph, &reference);

            //Write variations to file
            write_to_file(&PathBuf::from(out_path), &variations).unwrap();
        } else {
            panic!("Couldn't parse gfa file!");
        }
    } else {
        panic!("Couldn't parse reference file!");
    }
}

fn build_graph(gfa: &GFA) -> Graph<&str, &str> {
    let mut graph = Graph::<&str, &str>::new();

    let mut nodes_id = vec![];

    for seg in &gfa.segments {
        let id = graph.add_node(&seg.name);
        nodes_id.push(id);
    }

    let mut edges = vec![];
    for link in &gfa.links {
        let e = (&link.from_segment, &link.to_segment);
        edges.push(e);
    }

    println!("{:?}",edges);

    graph
}

fn find_variants(graph: &Graph<&str, &str>, reference: &String) -> Vec<Variant> {
    
    // TODO: function
    let mut variations: Vec<Variant> = Vec::new();
    
    //is this correct?
    //let first_node = graph.get_handle(NodeId::from(1), false);

    //traverse the graph
    //graph.for_each_handle(|h| {
    //   println!();
    //);

    variations
}

fn parse_reference(path: &PathBuf) -> Option<String> {
    let file = File::open(path).expect(&format!("Error opening file {:?}", path));
    
    let reader = BufReader::new(file);
    let mut reference = String::new();

    for line in reader.lines() {
        reference.push_str(&line.unwrap());
    }

    Some(reference)
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
