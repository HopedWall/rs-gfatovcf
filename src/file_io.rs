use chrono::Utc;
use json::*;
use std::fs::File;
use std::io::BufWriter;
use std::io::Write;
use std::path::PathBuf;

use crate::variant_identification::Variant;
use handlegraph::handlegraph::{handles_iter, HandleGraph};
use handlegraph::hashgraph::HashGraph;

/// Write variants to a VCF file
pub fn write_variants_to_file(
    path: &PathBuf,
    paths_list: &mut Vec<String>,
    variants: &[Variant],
) -> std::io::Result<()> {
    let mut file = File::create(path).unwrap_or_else(|_| panic!("Error creating file {:?}", path));

    // Add the .fa extension to paths if not present
    let mut paths_to_fasta: Vec<String> = paths_list
        .iter()
        .map(|x| {
            let mut var: String = x.clone();
            if !var.ends_with(".fa") && !var.ends_with(".fasta") {
                var.push_str(".fa");
            }
            format!("##reference={}", var)
        })
        .collect();
    // Create header
    let mut header: Vec<String> = Vec::new();
    header.push("##fileformat=VCFv4.2".to_string());
    header.push(format!(
        "##fileDate={}",
        Utc::now().format("%Y-%m-%d %H:%M:%S").to_string()
    ));
    header.append(&mut paths_to_fasta);
    header.push("##INFO=<ID=TYPE,Number=A,Type=String,Description=\"Type of each allele (snv, ins, del, mnp, complex)\">".to_string());
    header.push("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">".to_string());
    header.push(
        [
            "#CHROM",
            "POS",
            "ID",
            "REF",
            "ALT",
            "QUAL",
            "FILTER",
            "INFO",
            "FORMAT",
            "SampleName",
        ]
        .join("\t"),
    );

    file.write_all(header.join("\n").as_bytes())
        .expect("Error writing header");

    file.write_all(b"\n").expect("Error writing to file");
    for var in variants {
        file.write_all(var.to_vcf_string().as_bytes())
            .expect("Error writing variant");
    }
    Ok(())
}

/// Returns the json of a given HashGraph. A HashGraph is represented a JSON object
/// that contains a list of nodes and a list of edges (where each edge is represented
/// as a JSON Object with the start_node_id and end_node_id)
pub fn graph_to_json(graph: &HashGraph) -> JsonValue {
    //Create empty json
    let mut graph_json = JsonValue::new_object();

    // Obtains nodes and edges in the graph
    let mut nodes: Vec<u64> = handles_iter(graph).map(|x| u64::from(x.id())).collect();
    let mut edges: Vec<(u64, u64)> = std::iter::from_fn(graph.edges_iter_impl())
        .map(|edge| (u64::from(edge.0.id()), u64::from(edge.1.id())))
        .collect();

    // Sort both vecs so that they're easier to read
    nodes.sort();
    edges.sort();

    // Create an array of objects with from/where for each edge
    let mut edges_array = JsonValue::new_array();
    for (start, end) in edges {
        let temp = object! {
            start: start,
            end: end
        };
        edges_array.push(temp).unwrap();
    }

    graph_json["nodes"] = nodes.into();
    graph_json["edges"] = edges_array;

    graph_json
}

/// Writes a given JSON to a file
pub fn json_to_file(json: &JsonValue, path: &PathBuf) -> std::io::Result<()> {
    let mut buffer = BufWriter::new(File::create(path)?);
    json.write_pretty(&mut buffer, 4)?;
    Ok(())
}
