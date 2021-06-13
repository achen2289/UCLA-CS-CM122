from collections import defaultdict


class tree_node:
	def __init__(self, idx, color, children): # 0 = gray, 1 = red, 2 = blue, 3 = purple
		self.idx = idx
		self.color = color
		self.children = children


def form_tree(parent_to_children, node_colors):
	color_map = {"red": 1, "blue": 2}
	for node, color in node_colors.items():
		node_colors[node] = color_map[color]

	all_nodes = {}
	for node_idx, children in parent_to_children.items():
		curr_node = tree_node(node_idx, node_colors[node_idx], children)
		all_nodes[node_idx] = curr_node

	all_children = set()
	for node in all_nodes.values():
		for child in node.children:
			all_children.add(child)

	root = None
	for node_id in all_nodes.keys():
		if node_id not in all_children:
			root = all_nodes[node_id]
			break

	return root, all_nodes


# Color all children starting from given root
def color_tree(root, all_nodes): 
	if not root.children:
		return
	for child in root.children:
		color_tree(all_nodes[child], all_nodes)
	colors = list(set([all_nodes[child].color for child in root.children]))
	if len(colors) == 1:
		root.color = colors[0]
	else:
		root.color = 3


if __name__ == "__main__":
	with open("dataset_317419_6 (2).txt", "r") as f:
		content = f.read().strip("\n").split("\n")

	parent_to_children = {}
	node_colors = defaultdict(lambda: 0)
	reached_dash = False
	for line in content:
		if line == "-":
			reached_dash = True
			continue
		if not reached_dash:
			line = line.split(" -> ")
			if line[1] == "{}":
				children = []
			else:
				children = line[1].split(",")
			parent_to_children[line[0]] = children
		else:
			line = line.split(": ")
			node_colors[line[0]] = line[1]

	reverse_color_map = {0: "gray", 1: "red", 2: "blue", 3: "purple"}
	root, all_nodes = form_tree(parent_to_children, node_colors)
	color_tree(root, all_nodes)
	for i in range(len(all_nodes)):
		color = reverse_color_map[all_nodes[str(i)].color]
		print (f"{str(i)}: {color}")


