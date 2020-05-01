# BIOINF 7150 – Graph Genomes Prac

### Learning objectives

In this exercise you learn how to

- find toy examples to work with,
- construct a graph using `vg construct`,
- visualize that graph using `vg view`,
- simulate NGS reads from a graph using `vg sim`,
- create an index structure needed for read mapping using `vg index`,
- map reads to the graph using `vg map`.

### Getting started

First you will need to activate the conda environment that we will be using.

    conda activate variation

vg is a large collection of tools that is under very active development, to ensure we have access easily to the lastest versions of the software, we use the [Docker](https://www.docker.com/). It is already available on the course workstations. If you want to bulid it on your laptop, follow the instructions at the [vg homepage](https://github.com/vgteam/vg) (note that building vg and all submodules from source can take ~1h). 

To test if things are working smoothly, you can run:

    docker run quay.io/vgteam/vg date

The command above will run the vg "container" if it already exist on your machine, update it if there's a newer version, and download it if this is the first time you are to use the container.

In this exercise, you will use small toy examples from the `test` directory. So make sure you have checked out vg:

    cd data
    git clone https://github.com/vgteam/vg.git

Now create a directory to work on for this tutorial:

    mkdir exercise1
    cp vg/test/tiny/tiny* exercise1/
    cd exercise1

How do we give docker access to our files?

To simplify our lives, we are going to create an alias for docker's call to vg:

    alias vg="docker run -v $PWD:/data  quay.io/vgteam/vg vg"

The above command replaces the long command with `vg`. However, note that it binds always the local directory to `/data` path inside the container.

### Constructing and viewing your first graphs

Like many other toolkits, vg is comes with many different subcommands. First we will use `vg construct` to build our first graph. Run it without parameters to get information on its usage:

    vg construct

Let's construct a graph from just one sequence in file `tiny/tiny.fa`, which looks like this:

    >x
    CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG

To construct a graph, run

    vg construct -r tiny.fa -m 32 > tiny.ref.vg

This will create a (very boring) graph that just consists of a linear chain of nodes, each with 32 characters.

The switch `-m` tells vg to put at most 32 characters into each graph node. (You might want to run it with different values and observe the different results.) To visualize a graph, you can use `vg view`. Per default, `vg view` will output a graph in [GFA](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md) format. By adding `-j` or `-d`, you can generate [JSON](https://www.json.org/) or [DOT](https://www.graphviz.org/doc/info/lang.html) output.

    vg view /data/tiny.ref.vg
    vg view -d /data/tiny.ref.vg

To work with the JSON output the tool [jq](https://stedolan.github.io/jq/) comes in handy. To get all sequences in the graph, for instance, try

    # Install jq:
    conda install -c conda-forge jq

    vg view -j /data/tiny.ref.vg | jq '.node[].sequence'

Next, we use graphviz to layout the graph representation in DOT format.

    # Install dot:
    conda install -c anaconda graphviz

    vg view -d /data/tiny.ref.vg | dot -Tpdf -o tiny.ref.pdf

View the PDF and compare it to the input sequence. Now vary the parameter passed to `-m` of `vg construct` and visualize the result.

Ok, let's build a new graph that has some variants built into it. First, take a look at at `tiny/tiny.vcf.gz`, which contains variants in (gzipped) [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) format.

    vg construct -r /data/tiny.fa -v /data/tiny.vcf.gz -m 32 > tiny.vg

Visualize the outcome.  

Ok, that's nice, but you might wonder which sequence of nodes actually corresponds to the sequence (`tiny.fa`) you started from? To keep track of that, vg adds a **path** to the graph. Let's add this path to the visualization.

    vg view -dp /data/tiny.vg | dot -Tpdf -o tiny.pdf

You find the output too crowded? Option `-S` removes the sequence labels and only plots node IDs.

    vg view -dpS /data/tiny.vg | dot -Tpdf -o tiny.pdf

Another tool that comes with the graphviz package is *Neato*. It creates force-directed layouts of a graph.

    vg view -dpS /data/tiny.vg | neato -Tpdf -o tiny.pdf

For these small graphs, the difference is not that big, but for more involved cases, these layouts can be much easier to read.

A better solution for larger graphs is [Bandage](http://rrwick.github.io/Bandage/), which is designed for visualizing assembly graphs. It can't display path information, but it can scale to multi-megabase graphs with ease. Adjust the `--nodewidth` parameter (e.g. `--nodewidth 100`, up to 1000) if you are rendering a large graph and the resulting image appears blank.

    vg view /data/tiny.vg > x.gfa
    Bandage image /data/x.gfa /data/x.gfa.png

You can also apply `vg viz` to obtain a special kind of linear layout that scales well without losing path information.
It requires the `xg` index, which we'll build first.

    vg index -x /data/tiny.xg /data/tiny.vg
    vg viz -x /data/tiny.xg -o /data/tiny.svg

If the graph is very big, you'll need to view tiny.svg in chrome or chromium web browser.
