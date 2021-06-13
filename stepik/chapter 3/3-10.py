from typing import List, Dict, Tuple
from collections import defaultdict, namedtuple

Edge = namedtuple("Edge", ["start", "end"]) # represents an edge

def form_de_bruijn_graph(kmers: List[str]):
	graph = defaultdict(list) # map kmer string to node object
	prefix_to_node = defaultdict(list)

	for kmer in kmers:
		graph[kmer[:-1]].append(kmer[1:])

	return graph

def find_roots(graph):
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

def find_contigs_old(root, graph, all_contigs, curr_contig, visited_edges):
	children = graph[root]
	if len(children) == 0:
		all_contigs.append(curr_contig)
	elif len(children) == 1:
		if (root, children[0]) not in visited_edges:
			curr_contig += children[0][-1]
			visited_edges.append((root, children[0]))
			find_contigs(children[0], graph, all_contigs, curr_contig, visited_edges)
		else:
			all_contigs.append(curr_contig)
	elif len(children) > 1:
		all_contigs.append(curr_contig)
		for child in children:
			if (root, child) not in visited_edges:
				visited_edges.append((root, child))
				find_contigs(child, graph, all_contigs, root+child[-1], visited_edges)
			else:
				all_contigs.append(curr_contig)
		
def find_nonbranching_paths(graph, indegree_counts, outdegree_counts):
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
	all_contigs = []
	for path in all_paths:
		contig = path[0][0]
		for edge in path:
			contig += edge[1][-1]
		all_contigs.append(contig)
	return all_contigs

if __name__ == "__main__":
	file = "3-10-test.txt"
	file = "dataset_317292_5 (4).txt"
	# file = "Contigs/inputs/test3.txt"

	with open(file, "r") as f:
		kmers = f.read().strip("\n").split("\n")

	graph = form_de_bruijn_graph(kmers)
	roots, indegree_counts, outdegree_counts = find_roots(graph)

	print (f"indegree_counts: {indegree_counts}")
	print (f"outdegree_counts: {outdegree_counts}")

	# print (roots)
	# print (graph)

	# for root in roots:
	# 	contigs = find_contigs(root, graph, all_contigs, root, [])

	all_paths = find_nonbranching_paths(graph, indegree_counts, outdegree_counts)
	# print (all_paths)
	all_contigs = get_all_contigs(all_paths)
	
	print (" ".join(all_contigs))