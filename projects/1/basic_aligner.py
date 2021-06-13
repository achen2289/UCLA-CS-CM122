import sys
import argparse
import time
import zipfile
from collections import namedtuple, defaultdict

__author__ = "Alex Chen"

KMER_SIZE = 42
TOLERANCE = 5

SNP = namedtuple("SNP", ["Original", "SNP", "Position"])

def parse_reads_file(reads_fn):
    """
    :param reads_fn: the file containing all of the reads
    :return: outputs a list of all paired-end reads
    """
    try:
        with open(reads_fn, 'r') as rFile:
            print("Parsing Reads")
            first_line = True
            count = 0
            all_reads = []
            for line in rFile:
                count += 1
                if count % 1000 == 0:
                    print(count, " reads done")
                if first_line:
                    first_line = False
                    continue
                ends = line.strip().split(',')
                all_reads.append(ends)
        return all_reads
    except IOError:
        print("Could not read file: ", reads_fn)
        return None


def parse_ref_file(ref_fn):
    """
    :param ref_fn: the file containing the reference genome
    :return: a string containing the reference genome
    """
    try:
        with open(ref_fn, 'r') as gFile:
            print("Parsing Ref")
            first_line = True
            ref_genome = ''
            for line in gFile:
                if first_line:
                    first_line = False
                    continue
                ref_genome += line.strip()
        return ref_genome
    except IOError:
        print("Could not read file: ", ref_fn)
        return None

def index_genome(genome):
    index = defaultdict(list)
    size = KMER_SIZE // (TOLERANCE + 1)
    for i in range(len(genome)):
        if i + size <= len(genome):
            index[genome[i:i+size]].append(i)
    return index

def break_into_kmers(input_reads):
    kmers = []
    for read in input_reads:
        for i in range(len(read)):
            if i + KMER_SIZE <= len(read):
                kmers.append(read[i:i+KMER_SIZE])
    return kmers

def get_kmer_frequencies(kmers):
    frequencies = defaultdict(lambda: 0)
    for kmer in kmers:
        frequencies[kmer] += 1
    return frequencies

def filter_kmer_frequencies(frequencies):
    return {kmer: ct for kmer, ct in frequencies.items() if ct > 1}

def find_snps(kmers, genome_index, genome):
    def find_possible_starts():
        possible_starts = []
        section_sz = KMER_SIZE // (TOLERANCE + 1)
        # first, second, third = kmer[:section_sz], kmer[third_sz:2*section_sz], kmer[2*section_sz:]
        sections = []
        for i in range(TOLERANCE+1):
            sections.append(kmer[i*section_sz : (i+1)*section_sz])

        for i, section in enumerate(sections):
            starts = [start - i*section_sz for start in genome_index[section]]
            possible_starts.extend(starts)

        # possible_starts.extend(genome_index[first])

        # starts = [start - third_sz for start in genome_index[second]]
        # possible_starts.extend(starts)

        # starts = [start - 2*third_sz for start in genome_index[third]]
        # possible_starts.extend(starts)

        possible_starts = list(filter(lambda x: x>=0, possible_starts))

        return possible_starts

    def find_kmer_snps():
        for start in possible_starts:
            section = genome[start:start+KMER_SIZE]
            if len(section) != KMER_SIZE:
                continue
            curr_snps = defaultdict(list)
            for i in range(KMER_SIZE):
                if kmer[i] != section[i]:
                    snp = SNP(section[i], kmer[i], start+i)
                    curr_snps[start+i].append(snp)
                    if len(curr_snps) > 2:
                        break
            if len(curr_snps) <= 2:
                for location, present_snps in curr_snps.items():
                    all_snps[location].extend(present_snps)

    def concensus():
        new_all_snps = {}
        for location, location_snps in all_snps.items():
            max_snp, max_ct = None, 0
            counts = defaultdict(lambda: 0)
            for snp in location_snps:
                counts[snp[2]] += 1
                if counts[snp[2]] > max_ct:
                    max_snp = snp
                    max_ct = counts[snp[2]]
            if max_ct > 12:
                new_all_snps[max_snp[2]] = max_snp
        return new_all_snps

    all_snps = defaultdict(list)

    for kmer in kmers:
        possible_starts = find_possible_starts()
        possible_starts = list(set(possible_starts))
        find_kmer_snps()

    all_snps = concensus()

    return list(all_snps.values())

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='basic_aligner.py takes in data for homework assignment 1 consisting '
                                     'of a genome and a set of reads and aligns the reads to the reference genome, '
                                     'then calls SNPs based on this alignment')
    parser.add_argument('-g', '--referenceGenome', required=True, dest='reference_file',
                        help='File containing a reference genome.')
    parser.add_argument('-r', '--reads', required=True, dest='reads_file',
                        help='File containg sequencing reads.')
    parser.add_argument('-o', '--outputFile', required=True, dest='output_file',
                        help='Output file name.')
    parser.add_argument('-t', '--outputHeader', required=True, dest='output_header',
                        help='String that needs to be outputted on the first line of the output file so that the '
                             'online submission system recognizes which leaderboard this file should be submitted to.'
                             'This HAS to be practice_W_1_chr_1 for the practice data and hw1_W_2_chr_1 for the '
                             'for-credit assignment!')
    args = parser.parse_args()
    reference_fn = args.reference_file
    reads_fn = args.reads_file

    input_reads = parse_reads_file(reads_fn)
    if input_reads is None:
        sys.exit(1)
    reference = parse_ref_file(reference_fn)
    if reference is None:
        sys.exit(1)

    genome_index = index_genome(reference)
    # input_reads = [read for read_pair in input_reads for read in read_pair]
    new_input_reads = []
    new_input_reads.extend([read_pair[0] for read_pair in input_reads])
    new_input_reads.extend([read_pair[0][::-1] for read_pair in input_reads])
    new_input_reads.extend([read_pair[1] for read_pair in input_reads])
    new_input_reads.extend([read_pair[1][::-1] for read_pair in input_reads])
    input_reads = new_input_reads
    kmers = break_into_kmers(input_reads)
    kmer_frequencies = get_kmer_frequencies(kmers)
    kmer_frequencies = filter_kmer_frequencies(kmer_frequencies)
    kmers = list(kmer_frequencies.keys())
    snps = find_snps(kmers, genome_index, reference)

    # snps.sort(key=lambda snp: snp[2])

    output_fn = args.output_file
    zip_fn = output_fn + '.zip'
    with open(output_fn, 'w') as output_file:
        header = '>' + args.output_header + '\n>SNP\n'
        output_file.write(header)
        for x in snps:
            line = ','.join([str(u) for u in x]) + '\n'
            output_file.write(line)

        tails = ('>' + x for x in ('STR', 'CNV', 'ALU', 'INV', 'INS', 'DEL'))
        output_file.write('\n'.join(tails))

    with zipfile.ZipFile(zip_fn, 'w') as myzip:
        myzip.write(output_fn)
