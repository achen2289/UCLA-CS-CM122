from collections import defaultdict
from copy import deepcopy

def BWT(text):
	string = text
	rotations = []
	for i in range(len(text)):
		rotations.append(string)
		string = string[-1] + string[:-1]
	
	rotations = sorted(rotations)
	print (rotations)
	transform = [line[-1] for line in rotations]

	return "".join(transform)

def reconstruct(BWT):
	first_col = "".join(sorted(BWT))
	first_col_occur = defaultdict(dict) # map index to occurence
	last_col_index = defaultdict(dict) # map occurence to index

	for idx, ch in enumerate(first_col):
		first_col_occur[ch][idx] = len(first_col_occur[ch]) + 1

	for idx, ch in enumerate(BWT):
		last_col_index_size = len(last_col_index[ch])
		last_col_index[ch][last_col_index_size + 1] = idx

	# print (first_col_occur)
	# print (last_col_index)

	reconstructed = "$"
	curr_left_idx = last_col_index["$"][1]
	curr_left_ch = first_col[curr_left_idx]
	curr_occur = first_col_occur[curr_left_ch][curr_left_idx]
	curr_right_idx = last_col_index[curr_left_ch][curr_occur]
	curr_right_ch = BWT[curr_right_idx]

	while curr_right_idx != 0:
		reconstructed += curr_left_ch
		curr_left_idx = last_col_index[curr_left_ch][curr_occur]
		curr_left_ch = first_col[curr_left_idx]
		curr_occur = first_col_occur[curr_left_ch][curr_left_idx]
		curr_right_idx = last_col_index[curr_left_ch][curr_occur]
		curr_right_ch = BWT[curr_right_idx]
		
	reconstructed += curr_left_ch
	reconstructed = reconstructed[1:] + reconstructed[0]

	return reconstructed

def gen_last_to_first(first_col, BWT):
	occur = defaultdict(list)

	for i in reversed(range(len(first_col))):
		ch = first_col[i]
		occur[ch].append(i)

	last_to_first = defaultdict()
	for i, ch in enumerate(BWT):
		last_to_first[i] = occur[ch][-1]
		occur[ch].pop()

	return last_to_first


def BMWMatching(BWT, pattern, last_to_first):
	first_col = "".join(sorted(BWT))
	top, bottom = 0, len(BWT) - 1
	while top <= bottom:
		if pattern:
			sym = pattern[-1]
			pattern = pattern[:-1]
			top_idx, bottom_idx = None, None
			for i in range(top, bottom+1):
				if BWT[i] == sym:
					top_idx = i if not top_idx else top_idx
					bottom_idx = i
			if not top_idx or not bottom_idx:
				return 0
			top = last_to_first[top_idx]
			bottom = last_to_first[bottom_idx]
		else:
			return bottom - top + 1

def gen_first_occurence(BWT):
	first_occurence = defaultdict(lambda: 0)
	first_col = "".join(sorted(BWT))
	for i, ch in enumerate(first_col):
		if ch not in first_occurence:
			first_occurence[ch] = i
	return first_occurence

def get_counts(BWT):
	count = []
	count.append(defaultdict(lambda: 0))
	for i in range(len(BWT)):
		ch = BWT[i]
		ct = count[-1]
		count.append(deepcopy(ct))
		count[-1][ch] = count[-2][ch] + 1
	return count

def BetterBWMatching(first_occurence, BWT, pattern, count):
	top, bottom = 0, len(BWT) - 1
	while top <= bottom:
		if pattern:
			sym = pattern[-1]
			pattern = pattern[:-1]
			low_rank = count[top][sym]
			high_rank = count[bottom+1][sym]
			if high_rank > 0 and high_rank != low_rank:
				low_rank += 1
				top = first_occurence[sym] + low_rank - 1
				bottom = first_occurence[sym] + high_rank - 1
			else:
				return 0
		else:
			return bottom - top + 1

if __name__ == "__main__":
	# with open("dataset_317416_4 (4).txt", "r") as f:
	# 	content = f.read().strip("\n").split("\n")
	# 	orig = content[0]
	# 	patterns = content[1:]

	# # print (orig)
	# # print (patterns)
	# bwt = BWT(orig)
	# suffix_array = {}
	# for i in range(len(orig)):
	# 	suffix_array[orig[i:]] = i
	# first_occurence = gen_first_occurence(bwt)
	# count = get_counts(bwt)
	# results = []
	# for pattern in patterns:
	# 	if BetterBWMatching(first_occurence, bwt, pattern, count):
	# 		for suffix, idx in suffix_array.items():
	# 			if len(suffix) >= len(pattern) and suffix[:len(pattern)] == pattern:
	# 				results.append(idx)
	# results = sorted(results)
	# results = [str(result) for result in results]
	# print (" ".join(results))

	print (BWT("GATCCGC$"))
	print (reconstruct("CGGTC$CA"))

# if __name__ == "__main__":
# 	with open("dataset_317414_7 (1).txt", "r") as f:
# 		content = f.read().strip("\n").split("\n")
# 		bwt = content[0]
# 		patterns = content[1].split()

# 	# bwt = "GGCGCCGC$TAGTCACACACGCCGTA"
# 	# patterns = ["ACC", "CCG", "CAG"]

# 	# bwt = "hellomyname"
# 	first_occurence = gen_first_occurence(bwt)
# 	count = get_counts(bwt)


# 	results = []
# 	for pattern in patterns:
# 		results.append(str(BetterBWMatching(first_occurence, bwt, pattern, count)))

# 	print (" ".join(results))
# 	# print (BetterBWMatching(first_occurence, bwt, pattern, count))

# 	# count = get_counts("smnpbnnaaaaa$a")
# 	# for ct in count:
# 	# 	print (ct)


