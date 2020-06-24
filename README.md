# GFAtoVCF
GFAtoVCF is a tool that finds variants in a Variation Graph. Refer to [this link](https://gsocgraph.blogspot.com/) for further details.

## How to run this tool
First, clone this repo and move to the resulting directory:

```
git clone https://github.com/HopedWall/GFAtoVCF
cd GFAtoVCF
```

Create a new **dep** folder, that will contain the dependecies for this program. Then move to this folder.

```
mkdir dep
cd dep
```

Download the required dependencies via git clone, and go back to the main folder

```
git clone https://github.com/HopedWall/rs-gfa
git clone https://github.com/HopedWall/rs-handlegraph
cd ..
```

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

You can also add the **-v** (or **--verbose**) argument to display various debug messages.

## How it works
GFAtoVCF performs two main tasks during its execution:
1. **Bubble Detection**
2. **Variant Calling**

### Bubble Detection
A bubble consists of multiple directed unipaths (a path in which all internal vertices are of degree 2) between two vertices. Bubbles are commonly caused by a small number of errors in the centre of reads. More information can be found in [this paper](https://www.sciencedirect.com/science/article/pii/S0304397515009147#br0100)

Detecting bubbles is important for converting a file from GFA to VCF since Variants will be found exclusively inside bubbles. In order to detect bubbles, we fist convert the starting graph into a tree; then, by looking at how many nodes are at the each level of the tree, we can easily understand whether or not a bubble has been found. The approach will now be described in more detail:

1. **Run a BFS on the graph starting from an arbitrary node X**, by doing so we obtain a tree rooted in X. We chose a BFS instead of a DFS since we wanted the order of the edges to not matter, and BFS achieves exactly that.
2. **Detect, for each level of the tree, how many nodes are at that level**, this will provide crucial information for detecting bubbles in the next step.
3. **Detect bubbles** by analyzing the data structure obtained in the previous step; basically, we explore the tree level-by-level. The following situations are possible:
   - no bubble has been detected yet, a level with only 1 node is found -> **OPEN the bubble**
   - a bubble has already been detected, a level with multiple nodes is found -> **EXTEND the bubble**
   - a bubble has already been detected, a level with only 1 node is found -> **CLOSE the bubble**
    
Each bubble is represented internally as a tuple: **(NodeId_Bubble_Start, NodeId_Bubble_End)**, this will allow us to easily explore the nodes in the graph inside the bubble. Note that this script is currently unable to detect *nested* bubbles, also known as **superbubbles**. 

### Variant Calling
Variant Calling consists in finding **Variants**, which are differences in terms of sequences from a reference genome. In the graph, all genomes (including the reference ones) are described as **paths**. This second steps mostly invoves exploring different paths and comparing them to the references; instead of comparing all possible paths on the graph (which would be computationally expensive), we only focus on paths inside each bubble (since Variants can **only** be found there).

- For each reference path (from the GFA file) and for each bubble:
    - find all possible paths from NodeId_Bubble_Start to NodeId_Bubble_End. **Note that, by how bubbles are found, both NodeId_Bubble_Start and NodeId_Bubble_End will be present in all reference paths.**
    - follow both paths (reference and other path) node by node, at the same time:
        - if reference and path cross different nodes -> **SNV**
        - if reference reaches NodeId_Bubble_End before path -> **INS**
        - if path reaches NodeId_Bubble_End before reference -> **DEL**

