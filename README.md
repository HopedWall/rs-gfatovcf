# GFAtoVCF
GFAtoVCF is a tool that finds variants in a Variation Graph. Refer to [this link](https://gsocgraph.blogspot.com/) for further details.

## How to run this tool
First, clone this repo and move to the resulting directory:

```
git clone --recursive https://github.com/HopedWall/GFAtoVCF
cd GFAtoVCF
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

Detecting bubbles is important for converting a file from GFA to VCF since Variants will be found exclusively inside bubbles. In order to detect bubbles, we first build a support tree from the starting graph; then, by looking at how many nodes are at the each level of the tree, we can easily understand whether or not a bubble has been found. The approach will now be described in more detail:

1. **Run a BFS on the graph starting from a node X with no incoming edges**, by doing so we obtain a tree rooted in X. We chose a BFS instead of a DFS since we wanted the order of the edges to not matter, and BFS achieves exactly that.
2. **Detect, for each level of the tree, how many nodes are at that level**, this will provide crucial information for detecting bubbles in the next step.
3. **Detect bubbles** by analyzing the data structure obtained in the previous step; basically, we explore the tree level-by-level. The following situations are possible:
   - no bubble has been detected yet, a level with only 1 node is found -> **OPEN the bubble**
   - a bubble has already been detected, a level with multiple nodes is found -> **EXTEND the bubble**
   - a bubble has already been detected, a level with only 1 node is found -> **CLOSE the bubble**
    
Each bubble is represented internally as a tuple of NodeIDs: **(Bubble_Start, Bubble_End)**, this will allow us to easily explore the nodes in the graph inside the bubble. Note that this script is currently unable to detect *nested* bubbles, also known as **superbubbles**. 

### Variant Calling
Variant Calling consists in finding **Variants**, which are differences in terms of sequences from a reference genome. In the graph, all genomes, including the reference ones, are described as **paths**. This second steps mostly invoves exploring different paths and comparing them to the references; instead of comparing all possible paths in the graph, we only focus on paths inside each bubble (since Variants can **only** be found there). The approach will now be described in more detail:

1. **Consider a reference path**: remember that a path is a list of nodes, so ref = \[r1,r2...rn\]
2. **Compute all possible paths between Bubble_Start and Bubble_End**: each path will be represented as: path = \[p1,p2...pm\]
3. **Compare reference path to each possible path**. Consider a generic index i s.t. i <= min(n,m):
    - if ref\[i\] and path\[i\] are equal -> no variant found
    - if ref\[i+1\] and path\[i\] are equal -> **DEL**
    - if ref\[i\] and path\[i+1\] are equal -> **INS**
    - if ref\[i\] and path\[i\] are different -> **SNV**
    
This will be repeated for each reference path and for each bubble.

