# GFAtoVCF
GFAtoVCF is a tool that finds variants in a Variation Graph. Refer to [this link](https://gsocgraph.blogspot.com/) for further details.

## How to run this tool
To compile the program, run:

```
cargo build --release
```

The compiled executable can be found in **./target/release**. In order to run it, use the following command:

```
./gfa_to_vcf -i {input_file.gfa} -o {output_file.vcf}
```

Where:
- **input_file.gfa** is the gfa file that will be used to build the graph
- **output_file.vcf** is the vcf file that will be used to store the variants

## How it works
- Run a BFS on the graph from a node X, obtain a tree rooted in X
- Explore the tree level-by-level to obtain bubbles:
    - when at a given level there is only a single node: **OPEN/CLOSE the bubble**
    - when there are multiple nodes: **EXTEND the bubble**
    
    Each bubble is representing as a tuple: (NodeId_Bubble_Start, NodeId_Bubble_End)
- For each reference path (from the GFA file) and for each bubble:
    - find all possible paths from NodeId_Bubble_Start to NodeId_Bubble_End. **Note that, by how bubbles are found, both NodeId_Bubble_Start and NodeId_Bubble_End will be present in all reference paths.**
    - follow both paths (reference and other path) node by node, at the same time:
        - if reference and path cross different nodes -> **SNV**
        - if reference reaches NodeId_Bubble_End before path -> **INS**
        - if path reaches NodeId_Bubble_End before reference -> **DEL**

