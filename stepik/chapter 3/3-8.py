from typing import Dict, List, Set
from collections import namedtuple, defaultdict
from random import choice
from pprint import pprint

algorithms = __import__("3-3")

def eulerian_cycle(graph: Dict[str, List[str]]):
	Edge = namedtuple("Edge", ["start", "end"]) # represents an edge
	cycle = [] # represents edges traversed
	curr_node = list(graph.keys())[0] # arbitrary start node
	unfinished_nodes, unfinished_node_idxs = set(), {curr_node: 0}
	unfinished_nodes.add(curr_node)

	while len(unfinished_nodes) > 0:
		# pprint (f"cycle: {cycle}")
		# print (f"unfinished_nodes: {unfinished_node_idxs}")

		curr_node = unfinished_nodes.pop()
		unfinished_node_idx = unfinished_node_idxs.pop(curr_node)
		next_cycle = cycle[unfinished_node_idx:] + cycle[:unfinished_node_idx]
		# print (f"unfinished_nodes: {unfinished_node_idxs}")
		for node, idx in unfinished_node_idxs.items():
			if int(idx) < unfinished_node_idx:
				unfinished_node_idxs[node] += len(next_cycle) - unfinished_node_idx
			else:
				unfinished_node_idxs[node] -= unfinished_node_idx

		# print (f"unfinished_nodes: {unfinished_node_idxs}")

		# pprint (f"next_cycle: {next_cycle}")
		# print ("\n")

		while len(graph[curr_node]) > 0:
			next_node = graph[curr_node].pop()
			if len(graph[curr_node]) > 0:
				unfinished_nodes.add(curr_node)
				unfinished_node_idxs[curr_node] = len(next_cycle)
			else:
				if curr_node in unfinished_nodes:
					unfinished_nodes.discard(curr_node)
					unfinished_node_idxs.pop(curr_node)

			edge = Edge(curr_node, next_node)
			next_cycle.append(edge)
			curr_node = next_node

		cycle = next_cycle

	# pprint (f"cycle: {cycle}")
	# pprint (f"next_cycle: {next_cycle}")

	return cycle

def find_euler_cycle(graph):
	return eulerian_cycle(graph)

def find_euler_path(graph, added_edge):
	one, two = added_edge
	euler_cycle = eulerian_cycle(graph)
	for idx, edge in enumerate(euler_cycle):
		if edge[0] == one and edge[1] == two or edge[0] == two and edge[1] == one:
			euler_path = euler_cycle[idx+1:] + euler_cycle[:idx]
			break

	return euler_path

def form_graph_from_edges(file):
	graph = defaultdict(list)

	with open(file, "r") as f:
		content = f.read().strip("\n").split("\n")
		for line in content:
			line = line.split(" -> ")
			node = line[0]
			neighbors = line[1].split(",")
			for neighbor in neighbors:
				graph[node].append(neighbor)

	return graph

def form_graph_from_kmers(file, test_kmers=None):
	graph = defaultdict(list)

	with open(file, "r") as f:
		content = f.read().strip("\n").split("\n")
		kmers = content[1:]

	kmers = test_kmers

	graph = algorithms.form_debruijn_graph(kmers)
	return graph

def form_k_universal_string(file):
	with open(file, "r") as f:
		content = f.read().strip("\n")
		k = int(content)

	kmers = []
	for i in range(2 ** k):
		kmers.append(format(i, 'b').zfill(k))

	graph = algorithms.form_debruijn_graph(kmers)
	return k, graph

if __name__ == "__main__":
	# file = "StringReconstruction/inputs/test1.txt"
	file = "dataset_317290_11 (1).txt"
	# file = "k-universal-test.txt"

	kmers = []

	text = "ACTCACTCGACTT"
	k = 3

	for i in range(len(text)-k+1):
		kmers.append(text[i:i+k])

	# graph = form_graph_from_edges(file)
	graph = form_graph_from_kmers(file, kmers)
	# k, graph = form_k_universal_string(file)
	print (f"graph: {graph}")

	indegree_counts = defaultdict(lambda: 0)
	outdegree_counts = defaultdict(lambda: 0)
	for node, neighbors in graph.items():
		outdegree_counts[node] += len(neighbors)
		for neighbor in neighbors:
			indegree_counts[neighbor] += 1

	print (f"outdegree_counts: {outdegree_counts}")
	print (f"indegree_counts: {indegree_counts}")

	unbalanced_nodes = []
	for degree_counts in (indegree_counts, outdegree_counts):
		for node in degree_counts:
			if indegree_counts[node] != outdegree_counts[node]:
				unbalanced_nodes.append(node) if node not in unbalanced_nodes else None

	print (f"unbalanced_nodes: {unbalanced_nodes}")

	euler_path = None
	if len(unbalanced_nodes) == 2:
		one, two = unbalanced_nodes[0], unbalanced_nodes[1]
		if indegree_counts[one] < outdegree_counts[one]:
			graph[two].append(one)
		else:
			graph[one].append(two)
		euler_path = find_euler_path(graph, (one, two))
		# print (euler_path)
	else:
		euler_cycle = find_euler_cycle(graph)
	
	path = euler_path if euler_path else euler_cycle
	print (path)

	# reconstruction = algorithms.edges_to_genome(path)
	# print (reconstruction)

	# output = f"{path[0][0]}->"
	# output += "->".join([edge[1] for edge in path])

	# output = f"{path[0][0] + path[0][1][0]}"
	# output += "".join([edge[1][-1] for edge in path[:-(k)]])

	print (output)


