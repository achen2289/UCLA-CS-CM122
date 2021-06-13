from typing import List, Dict, Tuple
from collections import defaultdict

class kmer_node_old:
	def __init__(self, kmer: str):
		self.kmer = kmer
		self.children = []

class kmer_node:
	def __init__(self, kmer: str):
		self.label = kmer[:-1]
		self.kmer = kmer
		self.children = []

def form_graph_no_dup(kmers: List[str]) -> (kmer_node, Dict[str, kmer_node]):
	graph = {} # map kmer string to node object
	prefix_to_node = defaultdict(list)
	all_children = set()
	for kmer in kmers:
		graph[kmer] = kmer_node(kmer)
		prefix_to_node[kmer[:-1]].append(graph[kmer])

	for kmer, node in graph.items():
		if kmer[1:] in prefix_to_node:
			node.children = prefix_to_node[kmer[1:]]
			for child in node.children:
				all_children.add(child.kmer)

	for kmer in graph:
		if kmer not in all_children:
			root = graph[kmer]
			break

	return root, graph

def form_debruijn_graph(kmers: List[str]):
	graph = defaultdict(list) # map kmer string to node object
	prefix_to_node = defaultdict(list)
	all_children = set()
	for kmer in kmers:
		graph[kmer[:-1]].append(kmer[1:])

	return graph

def path_to_genome(path) -> str:
	return path[0] + "".join([comp[-1] for comp in path[1:]]) if path else None

def edges_to_genome(path: List[Tuple[str]]):
	return path[0][0] + "".join([comp[-1][-1] for comp in path]) if path else None

def four_universal_str():
	curr_num = 0
	curr_str = ""
	while curr_num < (2 ** 4):
		curr_str += format(curr_num, "b").zfill(4)
		curr_num += 1
	print (curr_str)

if __name__ == "__main__":
	test_file = "test_kmer.txt"
	actual_file = "dataset_317287_8.txt"

	file = actual_file

	# with open (file, "r") as f:
	# 	path = f.read().strip("\n").split("\n")
	# genome = path_to_genome(path)
	# print (genome)

	with open (file, "r") as f:
		content = f.read().strip("\n").split("\n")
		# k = int(content[0])
		# text = content[1]
		kmers = content

	kmers = []

	text = "ACTCACTCGACTT"
	k = 3

	for i in range(len(text)-k+1):
		kmers.append(text[i:i+k])

	graph = form_debruijn_graph(kmers)
	
	for kmer, neighbors in graph.items():
		print (f"{kmer} -> {','.join([neighbor for neighbor in neighbors])}")
	# for label, node in graph.items():
	# 	if node.children:
	# 		print (label, "->", ", ".join([child.label for child in node.children]))

	# four_universal_str()