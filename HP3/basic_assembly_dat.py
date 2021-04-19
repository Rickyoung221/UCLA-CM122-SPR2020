from os.path import join
import sys
import time
from collections import defaultdict, Counter
import sys
import os
import zipfile
import argparse
import random
sys.path.insert(0, os.path.abspath(".."))
sys.path.insert(0, os.path.abspath("../.."))


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


def preprocess_pairs(paired_reads, k, threshold):
    kmer_dict = defaultdict(int)
    read_length = len(paired_reads[0][0])
    kmers_per_read = read_length - k
    for pair in paired_reads:
        for i in range(kmers_per_read+1):
            kmer_dict[pair[0][i:i+k+1]] += 1
            kmer_dict[pair[1][i:i+k+1]] += 1
    clean_kmers = []
    for key, val in kmer_dict.items():
        # threshold = random.choice([1, 2])
        if val > threshold:
            clean_kmers.append(key)
    return(clean_kmers)


def contig(kmers):
    def inDegree(adj, v):
        return sum(1 if v == u else 0 for value in list(adj.values()) for u in value)

    def outDegree(adj, v):
        return len(adj[v])

    def nonBranchNode(adj, v):
        return (inDegree(adj, v) == 1) and (outDegree(adj, v) == 1)

    # db graph construction
    adj = defaultdict(list)
    for kmer in kmers:
        adj[kmer[:-1]].append(kmer[1:])

    # generate contigs
    contigs = []
    starts = [v for v in adj if not nonBranchNode(adj, v)]    # starts
    for start in starts:
        for v in adj[start]:
            nextV = v
            path = [start, nextV]  # take any path
            while nonBranchNode(adj, nextV):
                nextV = adj[nextV][0]
                path.append(nextV)
            r = path[0]
            r += ''.join(p[-1] for p in path[1:])
            contigs.append(r)
    # ' '.join(sorted(contigs))
    return sorted(contigs)


def genome_path(path):
    return ''.join([e[0] for e in path])+path[-1][1:]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='basic_assembly.py takes in data for homework assignment 3 consisting '
                                     'of a set of reads and aligns the reads to the reference genome.')
    parser.add_argument('-r', '--reads', required=True, dest='reads_file',
                        help='File containg sequencing reads.')
    parser.add_argument('-o', '--outputFile', required=True, dest='output_file',
                        help='Output file name.')
    parser.add_argument('-t', '--outputHeader', required=True, dest='output_header',
                        help='String that needs to be outputted on the first line of the output file so that the\n'
                             'online submission system recognizes which leaderboard this file should be submitted to.\n'
                             'This HAS to be one of the following:\n'
                             '1) spectrum_A_1_chr_1 for 10K spectrum reads practice data\n'
                             '2) practice_A_2_chr_1 for 10k normal reads practice data\n'
                             '3) hw3all_A_3_chr_1 for project 3 for-credit data\n')
    args = parser.parse_args()
    reads_fn = args.reads_file

    input_reads = parse_reads_file(reads_fn)
    if input_reads is None:
        sys.exit(1)

    k = 37
    frequency_threshold = 2

    clean_kmers = preprocess_pairs(input_reads, k, frequency_threshold)
    result = contig(clean_kmers)
    contigs = sorted(map(genome_path, result))

    output_fn = args.output_file
    zip_fn = output_fn + '.zip'
    with open(output_fn, 'w') as output_file:
        output_file.write('>' + args.output_header + '\n')
        output_file.write('>ASSEMBLY\n')
        output_file.write('\n'.join(contigs))
    with zipfile.ZipFile(zip_fn, 'w') as myzip:
        myzip.write(output_fn)
