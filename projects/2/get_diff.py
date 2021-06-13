import difflib as dl

def get_diff(seq, s1, s2):
	for tag, i1, i2, j1, j2 in seq.get_opcodes(): 
	    if tag == "insert": 
	    	print (f"{tag} at {i1}: {s2[j1:j2]}")
	    elif tag == "delete":
	    	print (f"{tag} at {i1}: {s1[i1:i2]}")

if __name__ == "__main__":
	s1 = "AGCT"
	s2 = "AGC"
	seq = dl.SequenceMatcher(None, s1, s2)
	get_diff(seq, s1, s2)