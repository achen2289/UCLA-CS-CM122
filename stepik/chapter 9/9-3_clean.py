from collections import defaultdict

class tree_node:
    def __init__(self):
        self.children = {}
        self.valid = False

def get_suffixes(text):
	suffixes = defaultdict(list)
	for i in range(len(text)):
		suffixes[text[i]].append(text[i:])
	return suffixes

def build_tree(suffixes):
	root = tree_node()
	root.valid = True
	curr_node = root
	for base in suffixes:
		substrs = suffixes[base]
		if len(substrs) == 1:
			curr_node.children[substrs] = tree_node()
		else:
			start, end = 0, 0
			max_len = len(substrs[0])
			for substr in substrs:
				max_len = max(max_len, len(substr))
			while start < max_len and end < max_len:
				chs = []
				for substr in substrs:
					if start < len(substr) and end < len(substr):
						chs.append(substr[start:end+1])
				if len(set(chs)) > 1:
					curr_node.children[chs[0][:-1]] = tree_node()
					curr_node = curr_node.children[chs[0][:-1]]
					start = end
				else:
					end += 1


if __name__ == "__main__":
	suffixes = get_suffixes("ATAAATG$")
	build_tree(suffixes)
