from collections import namedtuple, defaultdict
from typing import List, Tuple, Dict

KDmer = namedtuple("KDmer", ["Kmer1", "Kmer2"])
Prefix = namedtuple("Prefix", ["Kmer1", "Kmer2"])
Suffix = namedtuple("Suffix", ["Kmer1", "Kmer2"])

def gen_kdmers(text: str, k: int, d: int) -> List[Tuple[KDmer]]:
	all_kdmers = []
	for i in range(len(text)):
		if i+k <= len(text) and i+k+d+k <= len(text):
			kmer1 = text[i:i+k]
			kmer2 = text[i+k+d:i+k+d+k]
			kdmer = KDmer(kmer1, kmer2)
			all_kdmers.append(kdmer)
	sorted_kdmers = sorted(all_kdmers, key=lambda x: x[0]+x[1])
	return sorted_kdmers

def gen_paired_de_bruijn_graph(kdmers: List[str]):
	graph = defaultdict(list)
	for kdmer in kdmers:
		kmer1, kmer2 = kdmer[0], kdmer[1]
		prefix = Prefix(kmer1[:-1], kmer2[:-1])
		suffix = Suffix(kmer1[1:], kmer2[1:])
		graph[prefix].append(suffix)
	return graph

def eulerian_cycle(graph: Dict[Tuple[KDmer], List[Tuple[KDmer]]]):
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

def reconstruction(d, path):
	res = path[0][0][0]
	for edge in path:
		res += edge[1][0][-1]
	res += path[-d-1][0][1]
	for edge in path[-d-1:]:
		res += edge[1][1][-1]
	return res

if __name__ == "__main__":
	file = "3-9-test.txt"
	file = "dataset_317291_16 (1).txt"
	# file = "PairedStringReconstruction/inputs/test1.txt"
	# text = "TAATGCCATGGGATGTT"
	# k = 3
	# d = 2
	# sorted_kdmers = gen_kdmers(text, k, d)
	# res = ""
	# for kdmer in sorted_kdmers:
	# 	res += f"({kdmer[0]}|{kdmer[1]})"

	# print (res)

	with open(file, "r") as f:
		content = f.read().strip("\n").split("\n")
		first_line = content[0].split(" ")
		k, d = int(first_line[0]), int(first_line[1])
		kdmers_unprocessed = content[1:]

	kdmers = []
	for kdmer_line in kdmers_unprocessed:
		kdmer_pair = kdmer_line.split("|")
		kdmer = KDmer(kdmer_pair[0], kdmer_pair[1])
		kdmers.append(kdmer)

	graph = gen_paired_de_bruijn_graph(kdmers)

	# for node in graph:
	# 	print (node, graph[node])

	# print (f"graph: {graph}\n")

	indegree_counts = defaultdict(lambda: 0)
	outdegree_counts = defaultdict(lambda: 0)
	for node, neighbors in graph.items():
		outdegree_counts[node] += len(neighbors)
		for neighbor in neighbors:
			indegree_counts[neighbor] += 1

	# print (f"outdegree_counts: {outdegree_counts}\n")
	# print (f"indegree_counts: {indegree_counts}\n")

	unbalanced_nodes = []
	for degree_counts in (indegree_counts, outdegree_counts):
		for node in degree_counts:
			if indegree_counts[node] != outdegree_counts[node]:
				unbalanced_nodes.append(node) if node not in unbalanced_nodes else None

	# print (f"unbalanced_nodes: {unbalanced_nodes}\n")

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
	# for edge in path:
	# 	print (edge)
	print ("path:")
	for edge in path:
		print (edge)
	print ("\n")

	reconstructed = reconstruction(d, path)
	print (reconstructed)
	# print ("GTGGTCGTGAGATGTTGA")

