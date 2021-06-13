from os.path import join
import sys
import time
from collections import defaultdict, namedtuple, Counter
from typing import List, Dict
import sys
import os
import zipfile
import argparse
sys.path.insert(0, os.path.abspath(".."))
sys.path.insert(0, os.path.abspath("../.."))

__author__ = "Alex Chen"

KMER_SIZE = 35
Edge = namedtuple("Edge", ["start", "end"]) # represents an edge

def parse_reads_file(reads_fn):
    """
    :param reads_fn: the file containing all of the reads
    :return: outputs a list of all paired-end reads
    """
    try:
        with open(reads_fn, 'r') as rFile:
            print("Parsing Reads")
            first_line = True
            count = 0
            all_reads = []
            for line in rFile:
                count += 1
                # if count % 1000 == 0:
                #     print(count, " reads done")
                if first_line:
                    first_line = False
                    continue
                ends = line.strip().split(',')
                all_reads.append(ends)
        return all_reads
    except IOError:
        print("Could not read file: ", reads_fn)
        return None

def break_into_kmers(input_reads):
    """Break the input_reads into kmers of length KMER_SIZE.

    """
    kmers = []
    for read in input_reads:
        for i in range(len(read)):
            if i + KMER_SIZE <= len(read):
                kmers.append(read[i:i+KMER_SIZE])
    return kmers

def get_kmer_frequencies(kmers):
    """Get the number of occurences of each kmer.

    """
    frequencies = defaultdict(lambda: 0)
    for kmer in kmers:
        frequencies[kmer] += 1
    return frequencies

def filter_kmer_frequencies(frequencies):
    """Filter kmers by frequency.

    """
    return {kmer: ct for kmer, ct in frequencies.items() if ct > 3}

def form_debruijn_graph(kmers: List[str]):
    """Form debruijn graph using kmers.

    """
    graph = defaultdict(list)
    prefix_to_node = defaultdict(list)

    for kmer in kmers:
        graph[kmer[:-1]].append(kmer[1:])

    return graph

def find_roots(graph):
    """Find the roots of the graph, as well as indegree and outdegree counts.

    """
    indegree_counts = defaultdict(lambda: 0)
    outdegree_counts = defaultdict(lambda: 0)
    for node in graph:
        outdegree_counts[node] += len(graph[node])
        for child in graph[node]:
            indegree_counts[child] += 1

    roots = []
    for node in graph:
        if indegree_counts[node] == 0:
            roots.append(node)

    return roots, indegree_counts, outdegree_counts

def find_nonbranching_paths(graph, indegree_counts, outdegree_counts):
    """Find the non-branching paths, representing contigs.

    """
    paths = []
    one_to_one_nodes = set()
    for node in graph:
        if not (indegree_counts[node] == 1 and outdegree_counts[node] == 1):
            if outdegree_counts[node] > 0:
                for child in graph[node]:
                    edge = Edge(node, child)
                    nonbranching_paths = [edge]
                    while indegree_counts[child] == 1 and outdegree_counts[child] == 1:
                        edge = Edge(child, graph[child][0])
                        nonbranching_paths.append(edge)
                        child = graph[child][0]
                    paths.append(nonbranching_paths)
        else:
            one_to_one_nodes.add(node)

    isolated_cycles = []
    while len(one_to_one_nodes) > 0:
        isolated_cycle = []
        node = one_to_one_nodes.pop()
        while graph[node][0] in one_to_one_nodes:
            edge = Edge(node, graph[node][0])
            isolated_cycle.append(edge)
            node = graph[node][0]
            one_to_one_nodes.discard(node)
        if len(isolated_cycle) > 0 and graph[node][0] == isolated_cycle[0][0]:
            edge = Edge(node, graph[node][0])
            isolated_cycle.append(edge)
            isolated_cycles.append(isolated_cycle)

    paths.extend(isolated_cycles)
    return paths

def get_all_contigs(all_paths):
    """Form contigs from the paths found.

    """
    all_contigs = []
    for path in all_paths:
        contig = path[0][0]
        for edge in path:
            contig += edge[1][-1]
        all_contigs.append(contig)
    return all_contigs


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='basic_assembly.py takes in data for homework assignment 3 consisting '
                                                 'of a set of reads and aligns the reads to the reference genome.')
    parser.add_argument('-r', '--reads', required=True, dest='reads_file',
                        help='File containg sequencing reads.')
    parser.add_argument('-o', '--outputFile', required=True, dest='output_file',
                        help='Output file name.')
    parser.add_argument('-t', '--outputHeader', required=True, dest='output_header',
                        help='String that needs to be outputted on the first line of the output file so that the\n'
                             'online submission system recognizes which leaderboard this file should be submitted to.\n'
                             'This HAS to be one of the following:\n'
                             '1) spectrum_A_1_chr_1 for 10K spectrum reads practice data\n'
                             '2) practice_A_2_chr_1 for 10k normal reads practice data\n'
                             '3) hw3all_A_3_chr_1 for project 3 for-credit data\n')
    args = parser.parse_args()
    reads_fn = args.reads_file

    input_reads = parse_reads_file(reads_fn)
    if input_reads is None:
        sys.exit(1)
   
    new_input_reads = []
    new_input_reads.extend([read for read_pair in input_reads for read in read_pair])
    new_input_reads.extend([read[::-1] for read_pair in input_reads for read in read_pair])
    input_reads = new_input_reads

    kmers = break_into_kmers(input_reads)
    kmer_frequencies = get_kmer_frequencies(kmers)
    kmer_frequencies = filter_kmer_frequencies(kmer_frequencies)
    kmers = list(kmer_frequencies.keys())
    print (len(kmers))

    graph = form_debruijn_graph(kmers)
    roots, indegree_counts, outdegree_counts = find_roots(graph)

    all_paths = find_nonbranching_paths(graph, indegree_counts, outdegree_counts)
    contigs = get_all_contigs(all_paths)

    output_fn = args.output_file
    zip_fn = output_fn + '.zip'
    with open(output_fn, 'w') as output_file:
        output_file.write('>' + args.output_header + '\n')
        output_file.write('>ASSEMBLY\n')
        output_file.write('\n'.join(contigs))
    with zipfile.ZipFile(zip_fn, 'w') as myzip:
        myzip.write(output_fn)
