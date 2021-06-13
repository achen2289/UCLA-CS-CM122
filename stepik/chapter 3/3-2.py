def kmer(k, text):
	kmers = []
	for i in range(len(text)):
		if i+k <= len(text):
			kmers.append(text[i:i+k])
	return kmers

if __name__ == "__main__":
	file = "dataset_317284_3.txt"
	with open (file, "r") as f:
		content = f.read().strip("\n").split("\n")
		k = int(content[0])
		text = content[1]
	kmers = kmer(k, text)
	print ("\n".join(kmers))