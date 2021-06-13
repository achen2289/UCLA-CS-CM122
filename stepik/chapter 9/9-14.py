def kmer(genome, pattern):
	length = len(pattern)
	breakdown = [pattern[:length//3], pattern[length//3:2*length//3], pattern[2*length//3:]]
	for sec_num, sec in enumerate(breakdown):
		for i in range(len(genome)):
			if i+len(sec) < len(genome) and genome[i:i+len(sec)] == sec:
				mismatch = 0
				if sec_num == 0:
					for j in range(length//3, length):
						if pattern[j] != genome[i+j]:
							mismatch += 1
					if mismatch > 2:
						return -1
					else:
						return i
				elif sec_num == 2:
					for j in range(2*length//3):
						if pattern[j] != genome[i-(2*length//3)+j]:
							mismatch += 1
					if mismatch > 2:
						return -1
					else:
						return i - (2*length//3)
				elif sec_num == 1:
					for j in range(length//3):
						if pattern[j] != genome[i-(length//3)+j]:
							mismatch += 1
					for j in range(2*length//3, length):
						if pattern[j] != genome[i-(length//3)+j]:
							mismatch += 1
					if mismatch > 2:
						return -1
					else:
						return i - length//3

if __name__ == "__main__":
	genome = "abcdef"
	patterns = ["accdae"]
	with open("dataset_317417_10 (4).txt", "r") as f:
		content = f.read().strip("\n").split("\n")
		genome, patterns, allowed_mismatches = content[0], content[1], int(content[2])

	patterns = patterns.split(" ")

	indices = []
	if allowed_mismatches == 2:
		for pattern in patterns:
			res = kmer(genome, pattern)
			if res != -1:
				indices.append(res)

	indices = sorted(indices)
	indices = [str(index) for index in indices]
	print (" ".join(indices))
