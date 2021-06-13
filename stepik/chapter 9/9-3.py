class trie_node:
    def __init__(self, id):
        self.children = {}
        self.valid = False
        self.id = id
        self.label = None

class mod_trie_node:
    def __init__(self, idx):
        self.children = {}
        self.valid = False
        self.idx = idx

class Trie:
    def __init__(self):
        self.root = trie_node(0)
        self.node_count = 1

    def insert(self, word: str) -> None:
        node = self.root
        node.valid = True
        for l in word:
            if l not in node.children:
                node.children[l] = trie_node(self.node_count)
                node.valid = True
                self.node_count += 1
            node = node.children[l]

    def insert_word(self, word: str, node) -> None:
        if word not in node.children:
            node.children[word] = trie_node(self.node_count)
            node.valid = True
            self.node_count += 1

class ModTrie:
    def __init__(self):
        self.root = mod_trie_node(-1)
        self.root.valid = True

    def insert(self, word: str) -> None:
        for i in range(len(word)):
        	curr_node = self.root
        	for j in range(i, len(word)):
        		curr_sym = word[j]
        		if curr_sym in curr_node.children:
        			curr_node = curr_node.children[curr_sym]
        		else:
        			curr_node.children[curr_sym] = mod_trie_node(j)
        			curr_node = curr_node.children[curr_sym]
        	if len(curr_node.children) == 0:
        		curr_node.label = i

def print_trie_1(curr_node):
	if curr_node.valid:
		for key, child in curr_node.children.items():
			print (f"{curr_node.id}->{child.id}:{key}")
			print_trie_1(child)

def main(dataset):
	# with open(dataset, "r") as f:
	# 	patterns = f.read().strip("\n").split("\n")
	# 	trie = form_trie(patterns)
	# 	print_trie_1(trie.root)
	with open(dataset, "r") as f:
		patterns = f.read().strip("\n")
		trie = form_mod_trie(patterns)
		# print_trie_1(trie.root)
	# idxs = []
	# for idx in range(len(patterns[0])):
	# 	if scan_trie(trie, patterns[0][idx:]):
	# 		idxs.append(str(idx))
	# print (" ".join(idxs))
	# suffixes = []
	# with open(dataset, "r") as f:
	# 	text = f.read().strip("\n")
	# 	for ch in reversed(text):
	# 		suffixes.append(ch) if len(suffixes) == 0 else suffixes.append(ch + suffixes[-1])

def form_trie(patterns):
	trie = Trie()
	for pattern in patterns:
		trie.insert(pattern)
	return trie

def form_mod_trie(string):
	trie = ModTrie()
	trie.insert(string)
	return trie
	
def form_trie_1(patterns):
	trie = {}
	for pattern in patterns:
		curr_node = trie
		for ch in pattern:
			if ch in curr_node:
				curr_node = curr_node[ch]
			else:
				curr_node[ch] = {}
				curr_node = curr_node[ch]
	return trie

def form_trie_2(patterns):
	trie = ({}, 0)
	node_num = 1
	for pattern in patterns:
		curr_node = trie[0]
		for ch in pattern:
			if ch in curr_node:
				curr_node = curr_node[ch][0]
			else:
				curr_node[ch] = ({}, node_num)
				curr_node = curr_node[ch][0]
				node_num += 1
	return trie

def scan_trie(trie, sequence):
	for ch in sequence:
		if trie[0].get(ch, None):
			trie = trie[0][ch]
			if not trie[0]:
				break
		else:
			return False
	return True if not trie[0] else False



def print_trie_2(root, prev_num):
	if root:
		for key, next_node in root.items():
			next_root, next_node_num = next_node
			print (f"{prev_num}->{next_node_num}:{key}")
			print_trie(next_root, next_node_num)

def trie_to_tree(trie):
	if len(trie) == 0:
		return

	for k, v in trie.items():
		print (k, v)
		trie_to_tree(v)
		if len(v) == 1:
			v_key = list(v.keys())[0]
			v_val = v[v_key]
			del trie[k]
			trie[k+v_key] = v_val
			# trie[k+v_key] = trie.pop(v_key)

# def form_tree(text):
# 	trie = {}
# 	for i in range(len(text)):
# 		curr_node = trie
# 		for j in range(i, len(text)):
# 			curr_sym = text[j]
# 			if curr_sym in curr_node:
# 				curr_node = curr_node[curr_sym]
# 			else:
# 				curr_node[curr_sym] = {}
# 				curr_node = curr_node[(curr_sym, j)]
# 		if len(curr_node) == 0:


if __name__ == "__main__":
	act_set = "dataset_317406_4 (8).txt"
	test_set = "test_str"
	main(test_set)