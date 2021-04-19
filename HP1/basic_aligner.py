import sys
import argparse
import time
import zipfile
import textwrap
import random
from collections import defaultdict, Counter, OrderedDict


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


"""
    TODO: Use this space to implement any additional functions you might need

"""


def align(reads, reference, k):
    """
    input:  array of reads, reference genome
    output: list of snps
    """
    def check_snp(read, genome_seq):
        index = None
        mismatches = 0
        for i in range(len(read)):
            if read[i] != genome_seq[i]:
                mismatches += 1
                index = i
        if 3 > mismatches:
            return[genome_seq[index], read[index], index]
        return None

    def get_snps():
        snps = []
        read_length = len(reads[0])
        for read in reads:
            mismatches = 0
            kmers = textwrap.wrap(read, k)
            first_match = None
            first_match_index = None
            for i, kmer in enumerate(kmers):
                if genome.get(kmer):
                    if not first_match:
                        first_match = kmer
                        first_match_index = i*k
                else:
                    mismatches += 1
            if mismatches > 0 and 2 > mismatches:
                index_array = genome.get(first_match)
                for i in index_array:
                    start = i - first_match_index
                    genome_seq = reference[start: start+read_length]
                    x = check_snp(read, genome_seq)
                    if x is not None:
                        snps.append([x[0], x[1], start+x[2]])
        return snps

    def create_ref_dict(k):
        d = OrderedDict()
        for i in range(len(reference)-k):
            seq = reference[i:i+k]
            d.setdefault(seq, []).append(i)
        return d

    genome = create_ref_dict(k)
    snps = get_snps()
    snps.sort(key=lambda x: x[2])

    # count only the snps that occur above a specific frequency
    snps_new = []
    cur_snp = snps[0]
    cur_snp_count = 1
    for s in snps[1:]:
        if str(s) == str(cur_snp):
            cur_snp_count += 1
        else:
            if cur_snp_count > 2:
                snps_new.append(cur_snp)
            cur_snp = s
            cur_snp_count = 1

    print(len(snps_new))
    return(snps_new)


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

    # paired reads into list of reads
    reads = []
    for pair in input_reads:
        reads.append(pair[0])
        reads.append(pair[1])

    k = 25
    snps = align(reads, reference, k)

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
