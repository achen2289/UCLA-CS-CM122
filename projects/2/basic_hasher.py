import sys
import argparse
import time
import zipfile
from collections import namedtuple, defaultdict
import difflib as dl

"""
Finds SNPS and INDELS between a genome and paired reads (FASTA format)

Read and genome parsing was provided: https://github.com/rosie068/CM122_starter_code
"""

__author__ = "Alex Chen"

READ_LENGTH = 50

KMER_SIZE = 38 # Used for finding SNPs
TOLERANCE = 2 # Error tolerance for SNPs
SNP_CONCENSUS_FACTOR = 18

INDEL_MAX_SIZE = 3 # Max number of BPs inserted or deleted from a section
KMER_SECTIONS = 5 # Used for finding INDELs
INDEL_CONCENSUS_FACTOR = 4 

SNP = namedtuple("SNP", ["Original", "SNP", "Position"])

def parse_reads_file(reads_fn):
    """
    :param reads_fn: the file containing all of the reads
    :return: outputs a list of all paired-end reads
    """
    try:
        with open(reads_fn, 'r') as rFile:
            print ("Parsing Reads")
            first_line = True
            count = 0
            all_reads = []
            for line in rFile:
                count += 1
                # if count % 1000 == 0:
                    # print(count, " reads done")
                if first_line:
                    first_line = False
                    continue
                ends = line.strip().split(',')
                all_reads.append(ends)
        return all_reads
    except IOError:
        print ("Could not read file: ", reads_fn)
        return None

def parse_ref_file(ref_fn):
    """
    :param ref_fn: the file containing the reference genome
    :return: a string containing the reference genome
    """
    try:
        with open(ref_fn, 'r') as gFile:
            print ("Parsing Ref")
            first_line = True
            ref_genome = ''
            for line in gFile:
                if first_line:
                    first_line = False
                    continue
                ref_genome += line.strip()
        return ref_genome
    except IOError:
        print ("Could not read file: ", ref_fn)
        return None

def index_genome(genome):
    """Index the genome, by sequences of length size_for_kmers and size_for_indels.

    """
    index_for_kmers, index_for_indels = defaultdict(list), defaultdict(list)
    size_for_kmers = KMER_SIZE // (TOLERANCE + 1)
    size_for_indels = READ_LENGTH // KMER_SECTIONS
    for i in range(len(genome)):
        if i + size_for_kmers <= len(genome):
            index_for_kmers[genome[i:i+size_for_kmers]].append(i)
        if i + size_for_indels <= len(genome):
            index_for_indels[genome[i:i+size_for_indels]].append(i)
    return index_for_kmers, index_for_indels

def break_into_kmers(input_reads):
    """Break the input_reads into kmers of length KMER_SIZE.

    """
    kmers = []
    for read in input_reads:
        for i in range(len(read)):
            if i + KMER_SIZE <= len(read):
                kmers.append(read[i:i+KMER_SIZE])
    return kmers

def get_kmer_frequencies(kmers):
    """Get the number of occurences of each kmer.

    """
    frequencies = defaultdict(lambda: 0)
    for kmer in kmers:
        frequencies[kmer] += 1
    return frequencies

def filter_kmer_frequencies(frequencies):
    """Filter kmers by frequency.

    """
    return {kmer: ct for kmer, ct in frequencies.items() if ct > 2}

def find_snps(kmers, genome_index, genome):
    """Given all kmers, find the SNPS.

    Break kmers up into TOLERANCE+1 pieces, then use index to find their position in the genome.

    """
    def find_possible_starts():
        possible_starts = []
        section_sz = KMER_SIZE // (TOLERANCE + 1)
        sections = []
        for i in range(TOLERANCE+1):
            sections.append(kmer[i*section_sz : (i+1)*section_sz])

        for i, section in enumerate(sections):
            starts = [start - i*section_sz for start in genome_index[section]]
            possible_starts.extend(starts)

        possible_starts = list(filter(lambda x: x>=0, possible_starts))

        return possible_starts

    def align_reads():
        """Narrow down kmers to those that have TOLERANCE or less mismatches.

        """
        aligned_read_locations = set()
        for start in possible_starts:
            # aligned_read_location = genome[start:start+KMER_SIZE]
            if start+KMER_SIZE > len(genome):
                continue
            # if len(aligned_read_location) < KMER_SIZE:
            #     continue
            total_mismatches = sum(kmer[i] != genome[start+i] for i in range(KMER_SIZE))
            if total_mismatches > 0 and total_mismatches <= TOLERANCE:
                aligned_read_locations.add(start)
        return aligned_read_locations

    def find_kmer_snps():
        """Use aligned read locations to determine the actual SNPS.

        """
        for start in aligned_read_locations:
            # section = genome[start:start+KMER_SIZE]
            # if len(section) != KMER_SIZE:
            #     continue
            # curr_snps = defaultdict(list)
            for i in range(KMER_SIZE):
                if kmer[i] != genome[start+i]:
                    snp = SNP(genome[start+i], kmer[i], start+i)
                    all_snps[start+i].append(snp)
                    # curr_max_snp, curr_max_ct = max_snps[start+i]
                    # if curr_max_snp and curr 
                    # if len(curr_snps) > 2:
                    #     break
            # if len(curr_snps) <= 2:
            #     for location, present_snps in curr_snps.items():
            #         all_snps[location].extend(present_snps)

    def concensus():
        """Run a concensus algorithm, keeping only the SNPs for each location 
        which have at least SNP_CONCENSUS_FACTOR kmers agreeing with that SNP.

        """
        new_all_snps = {}
        for location, location_snps in all_snps.items():
            max_snp, max_ct = None, 0
            counts = defaultdict(lambda: 0)
            for snp in location_snps:
                counts[snp[1]] += 1
                if counts[snp[1]] > max_ct:
                    max_snp = snp
                    max_ct = counts[snp[1]]
            if max_ct > SNP_CONCENSUS_FACTOR:
                new_all_snps[max_snp[2]] = max_snp
        return new_all_snps

    all_snps = defaultdict(list)

    for kmer in kmers:
        possible_starts = find_possible_starts()
        aligned_read_locations = align_reads()
        find_kmer_snps()
    
    all_snps = concensus()

    return list(all_snps.values()), all_snps

def genome_mod(reference, snps):
    """Update the reference genome with the determine SNPs to increase accuracy.

    """
    reference_arr = list(reference)
    for snp in snps:
        original, snp, position = snp
        reference_arr[position] = snp
    reference = "".join(reference_arr)
    return reference

def match_outer_sections(genome_index, genome, input_reads):
    """Match outer sections of each input read to genome, and return the 
    middle sections that have a valid length.

    """
    potential_indel_locations = []
    read_section_length = len(input_reads[0]) // KMER_SECTIONS
    middle_sec_base_size = int(((KMER_SECTIONS-1) / KMER_SECTIONS) * len(input_reads[0]))
    middle_sec_sizes = (middle_sec_base_size - INDEL_MAX_SIZE, middle_sec_base_size + INDEL_MAX_SIZE)

    for read in input_reads:
        first_sec = read[:read_section_length]
        last_sec = read[len(read) - read_section_length:]
        middle_sec = read[read_section_length:len(read)-read_section_length]
        if first_sec in genome_index and last_sec in genome_index:
            first_sec_locs = genome_index[first_sec]
            last_sec_locs = genome_index[last_sec]
            for last_loc in last_sec_locs:
                for first_loc in first_sec_locs:
                    if first_loc >= last_loc:
                        break
                    start_diff = last_loc - first_loc
                    if start_diff <= middle_sec_sizes[1] and start_diff >= middle_sec_sizes[0]:
                        start = first_loc + read_section_length
                        end = last_loc
                        potential_indel_locations.append((start, end, middle_sec))

    return potential_indel_locations

def find_indels(locations, genome):
    """Use the potential indel locations and run either smith waterman or use the 
    difflib to find INDELs.

    """
    insertions, deletions = defaultdict(list), defaultdict(list)
    for location in locations:
        genome_start, genome_end, read_section = location
        genome_section = genome[genome_start : genome_end]
        # smith_waterman(genome_section, read_section, genome_start, insertions, deletions)
        use_diff_lib(genome_section, read_section, genome_start, insertions, deletions)
    return insertions, deletions

def use_diff_lib(genome_section, read_section, genome_start, insertions, deletions):
    """Run difflib to get INDELs.

    """
    s1 = genome_section
    s2 = read_section
    seq_matcher = dl.SequenceMatcher(None, s1, s2)
    for tag, i1, i2, j1, j2 in seq_matcher.get_opcodes(): 
        if tag == "insert": 
            # print (f"{tag} at {i1}: {s2[j1 : j2]}")
            insertions[i1+genome_start].append((s2[j1 : j2], i1+genome_start))
        if tag == "delete":
            # print (f"{tag} at {i1}: {s1[i1:i2]}")
            deletions[i1+genome_start].append((s1[i1 : i2], i1+genome_start))

def smith_waterman(genome_section, read_section, genome_start, insertions, deletions):
    """Run smith waterman algorithm to get INDELs.

    """
    overall_max_score, overall_max_indices = 0, (0, 0)
    opt = [[0] * (len(read_section)+1) for i in range(len(genome_section)+1)]
    for i in range(1, len(genome_section)+1):
        curr_genome_char = genome_section[i-1]
        for j in range(1, len(read_section)+1):
            curr_read_char = read_section[j-1]
            max_score = max(opt[i-1][j] - 2, opt[i][j-1] - 2, opt[i-1][j-1] - 3)
            max_score = max(max_score, 0)
            if curr_genome_char == curr_read_char:
                max_score = max(max_score, opt[i-1][j-1] + 3)
            if max_score >= overall_max_score:
                overall_max_score = max_score
                overall_max_indices = (i, j)
    
            opt[i][j] = max_score

    r, c = overall_max_indices
    curr_ins, curr_del = "", ""
    while r > 0 or c > 0:
        if r-1 >= 0 and c-1 >= 0:
            diag, left, up = opt[r-1][c-1], opt[r][c-1], opt[r-1][c]
            if diag >= left and diag >= up:
                if curr_ins:
                    insertions.add((curr_ins[::-1], genome_start + r))
                    curr_ins = ""
                if curr_del:
                    deletions.add((curr_del[::-1], genome_start + r))
                    curr_del = ""
                r -= 1
                c -= 1
            elif left >= diag and left >= up:
                curr_ins += read_section[c-1]
                c -= 1
            elif up >= diag and up >= left:
                curr_del += genome_section[r-1]
                r -= 1
        elif c-1 >= 0:
            curr_ins += read_section[c-1]
            c -= 1
        elif r-1 >= 0:
            curr_del += genome_section[r-1]
            r -= 1

    if curr_ins:
        insertions.add((curr_ins[::-1], genome_start + r))
        curr_ins = ""
    if curr_del:
        deletions.add((curr_del[::-1], genome_start + r))
        curr_del = ""

def indel_concensus(insertions, deletions):
    """Run concensus algorithm in INDELs found.

    """
    new_insertions, new_deletions = {}, {}
    for loc, ins_at_loc in insertions.items():
        max_ins, max_ct = None, 0
        counts = defaultdict(lambda: 0)
        for ins in ins_at_loc:
            counts[ins[0]] += 1
            if counts[ins[0]] > max_ct:
                max_ins = ins
                max_ct = counts[ins[0]]
        if max_ct > INDEL_CONCENSUS_FACTOR:
            new_insertions[loc] = max_ins

    for loc, del_at_loc in deletions.items():
        max_del, max_ct = None, 0
        counts = defaultdict(lambda: 0)
        for dels in del_at_loc:
            counts[dels[0]] += 1
            if counts[dels[0]] > max_ct:
                max_del = dels
                max_ct = counts[dels[0]]
        if max_ct > INDEL_CONCENSUS_FACTOR:
            new_deletions[loc] = max_del

    insertions = list(new_insertions.values())
    deletions = list(new_deletions.values())

    return insertions, deletions

# Test smith waterman algorithm
if __name__ != "__main__":
    insertions, deletions = set(), set()
    smith_waterman("AGCTA", "AAGGCTA", 0, insertions, deletions)
    print (insertions, deletions)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='basic_hasher.py takes in data for homework assignment 2 consisting '
                                     'of a genome and a set of reads and aligns the reads to the reference genome, '
                                     'then calls SNPS and indels based on this alignment.')
    parser.add_argument('-g', '--referenceGenome', required=True, dest='reference_file',
                        help='File containing a reference genome.')
    parser.add_argument('-r', '--reads', required=True, dest='reads_file',
                        help='File containg sequencing reads.')
    parser.add_argument('-o', '--outputFile', required=True, dest='output_file',
                        help='Output file name.')
    parser.add_argument('-t', '--outputHeader', required=True, dest='output_header',
                        help='String that needs to be outputted on the first line of the output file so that the\n'
                             'online submission system recognizes which leaderboard this file should be submitted to.\n'
                             'This HAS to be one of the following:\n'
                             '1) practice_W_3_chr_1 for 10K length genome practice data\n'
                             '2) practice_E_1_chr_1 for 1 million length genome practice data\n'
                             '3) hw2undergrad_E_2_chr_1 for project 2 undergrad for-credit data\n'
                             '4) hw2grad_M_1_chr_1 for project 2 grad for-credit data\n')
    args = parser.parse_args()
    reference_fn = args.reference_file
    reads_fn = args.reads_file

    input_reads = parse_reads_file(reads_fn)
    if input_reads is None:
        sys.exit(1)
    reference = parse_ref_file(reference_fn)
    if reference is None:
        sys.exit(1)

    genome_index_for_snps, genome_index_for_indels = index_genome(reference)

    new_input_reads = []
    new_input_reads.extend([read for read_pair in input_reads for read in read_pair])
    new_input_reads.extend([read[::-1] for read_pair in input_reads for read in read_pair])
    input_reads = new_input_reads

    kmers = break_into_kmers(input_reads)
    kmer_frequencies = get_kmer_frequencies(kmers)
    kmer_frequencies = filter_kmer_frequencies(kmer_frequencies)
    kmers = list(kmer_frequencies.keys())
    snps, location_to_snps = find_snps(kmers, genome_index_for_snps, reference)

    print (f"SNP CT: {len(snps)}")
    print (f"FINISHED SNPS")

    reference = genome_mod(reference, snps)

    potential_indel_locations = match_outer_sections(genome_index_for_indels, reference, input_reads)
    insertions, deletions = find_indels(potential_indel_locations, reference)
    insertions, deletions = indel_concensus(insertions, deletions)

    print (f"INS CT: {len(insertions)}")
    print (f"DEL CT: {len(deletions)}")
    print (f"FINISHED INDELS")

    output_fn = args.output_file
    zip_fn = output_fn + '.zip'
    with open(output_fn, 'w') as output_file:
        output_file.write('>' + args.output_header + '\n>SNP\n')
        for x in snps:
            output_file.write(','.join([str(u) for u in x]) + '\n')
        output_file.write('>INS\n')
        for x in insertions:
            output_file.write(','.join([str(u) for u in x]) + '\n')
        output_file.write('>DEL\n')
        for x in deletions:
            output_file.write(','.join([str(u) for u in x]) + '\n')
    with zipfile.ZipFile(zip_fn, 'w') as myzip:
        myzip.write(output_fn)
