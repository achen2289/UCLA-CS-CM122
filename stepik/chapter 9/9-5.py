from collections import defaultdict

def main(input):
	all_words = defaultdict(lambda: 0)
	for i in range(len(input)):
		for j in range(i+1, len(input)):
			all_words[input[i:j+1]] += 1

	for word in sorted(all_words.keys(), key=lambda x: len(x), reverse=True):
		if all_words[word] > 1:
			return word 
	return None

def longest_common_substr(one, two):
	substrs_1, substrs_2 = defaultdict(lambda: 0), defaultdict(lambda: 0)
	for (text, substrs) in [(one, substrs_1), (two, substrs_2)]:
		all_words = substrs
		word = text
		for i in range(len(word)):
			for j in range(i+1, len(word)):
				all_words[word[i:j+1]] += 1

	print ("hi")

	sorted_1 = sorted(substrs_1.keys(), key=lambda x: len(x), reverse=True)
	sorted_2 = sorted(substrs_2.keys(), key=lambda x: len(x), reverse=True)
	print ("bye")
	for key in sorted_1:
		if key in substrs_2:
			return key
	return None

def shortest_common_substr(one, two):
	substrs_1, substrs_2 = defaultdict(lambda: 0), defaultdict(lambda: 0)
	for (text, substrs) in [(one, substrs_1), (two, substrs_2)]:
		all_words = substrs
		word = text
		for i in range(len(word)):
			for j in range(i+1, len(word)):
				all_words[word[i:j+1]] += 1

	print ("hi")

	all_sorted = sorted(substrs_1.keys(), key=lambda x: len(x))
	print ("bye")
	for key in all_sorted:
		if key not in substrs_2:
			return key
	return None

def longest_common_subseq(one, two):
	rows, cols = len(one) + 1, len(two) + 1
	opt = [[0] * cols] * rows
	for i in range(len(opt[0])):
		opt[0][i] = 0
	for i in range(len(opt)):
		opt[i][0] = 0
	for i in range(1, len(one)+1):
		for j in range(1, len(two)+1):
			if one[i] == two[j]:
				opt[i][j] = opt[i-1][j-1] + 1
			else:
				opt[i][j] = max(opt[i][j-1], opt[i-1][j])
	return opt[-1][-1]


if __name__ == "__main__":
	with open("dataset_317408_7 (2).txt", "r") as f:
		words = f.read().strip("\n").split("\n")
		print (shortest_common_substr(words[0], words[1]))
	# print (main("ATATCGTTTTATCGTT"))
