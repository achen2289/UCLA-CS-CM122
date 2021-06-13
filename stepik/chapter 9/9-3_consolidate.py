class trie_node:
    def __init__(self, id, edge=None):
        self.children = []
        self.edge = edge
        self.valid = False
        self.id = id

class Trie:
    def __init__(self):
        self.root = trie_node(0)
        self.node_count = 1

    def insert(self, word: str) -> None:
        node = self.root
        for l in word:
            next_node = None
            for child in node.children:
                if child.edge == l:
                    next_node = child
                    break
            if not next_node:
                next_node = trie_node(self.node_count, l)
                self.node_count += 1
                next_node.valid = True
                node.children.append(next_node)
            node = next_node

def trie_to_tree(node):
    if node.edge == "$":
        return
    for child in node.children:
        trie_to_tree(child)
        if len(child.children) == 1:
            child.edge += child.children[0].edge
            child.children = child.children[0].children

def build_trie(suffixes):
    trie = Trie()
    for suffix in suffixes:
        trie.insert(suffix)
    return trie

def print_trie_with_id(curr_node):
    for child in curr_node.children:
        print (f"{curr_node.id}->{child.id}:{child.edge}")
        print_trie_with_id(child)

def print_trie(node):
    if node:
        if node.valid:
            print (node.edge)
        for child in node.children:
            print_trie(child)

def get_suffixes(text):
    suffixes = []
    for i in range(len(text)):
        suffixes.append(text[i:])
    return suffixes

if __name__ == "__main__":
    suffixes = get_suffixes("ATAAATG$")
    # with open("dataset_317406_4 (9).txt", "r") as f:
    #     patterns = f.read().strip("\n").split("\n")
    with open("dataset_317408_4.txt", "r") as f:
        text = f.read().strip("\n")
    suffixes = get_suffixes(text)
    trie = build_trie(suffixes)
    # # print_trie_with_id(trie.root)
    trie_to_tree(trie.root)
    print_trie(trie.root)

