from os.path import join
import sys
import time
from collections import defaultdict, Counter, OrderedDict
import sys
import os
import zipfile
import argparse
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
    """
    :param  paired_reads: list of pairs of reads provided by input file
            k: size of kmer
            threshold: count as error and delete kmer if it does not occur more times than threshold  
    :return: list of kmers that occur above threshold frequency
    """
    kmers = []
    read_length = len(paired_reads[0][0])
    kmers_per_read = read_length - k
    for pair in paired_reads:
        for i in range(kmers_per_read+1):
            kmers.append(pair[0][i:i+k])
            kmers.append(pair[1][i:i+k])
    kmer_dict = Counter(kmers)
    clean_kmers = []
    for key, val in kmer_dict.items():
        if val > threshold:
            clean_kmers.append(key)
    return(clean_kmers)


def kmer_to_dbg(kmers):
    """
    :param  kmers: list of non-duplicated kmers
    :return: dictionary that is adjacency list representing de Bruijn graph
    """
    prefix_dict = OrderedDict()
    for kmer in kmers:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        prefix_dict.setdefault(prefix, []).append(suffix)
    for key, value in prefix_dict.items():
        prefix_dict[key] = sorted(value)
    return prefix_dict


def dbg_to_contigs(dbg):
    """
    :param  dbg: dictionary that is adjacency list representing de Bruijn graph
    :return: sorted list of contigs
    """
    def check_cycle(vertex):
        cycle = [vertex]
        while True:
            last = cycle[-1]
            if indegrees.get(last) != 1 or outdegrees.get(last) != 1:
                return None
            cycle.append(dbg.get(last)[0])
            if cycle[0] == last:
                return cycle
        return None

    paths = []
    indegrees = {}
    outdegrees = {}
    for key, value in dbg.items():
        outdegrees[key] = len(value)
        for v in value:
            indegrees[v] = indegrees.get(v, 0) + 1

    for key, value in dbg.items():
        indegree = indegrees.get(key)
        outdegree = outdegrees.get(key)
        if indegree == outdegree == 1:
            visited = False
            for path in paths:
                if key in path:
                    visted = True
                    break
            if not visited:
                cycle = check_cycle(key)
                if cycle:
                    paths.append(cycle)
        elif (outdegree > 0):
            for end_vertex in value:
                branch = [key, end_vertex]
                while indegrees.get(branch[-1]) == outdegrees.get(branch[-1]) == 1:
                    branch.append(dbg.get(branch[-1])[0])
                paths.append(branch)
    return sorted(map(lambda p: ''.join([i[0] for i in p])+p[-1][1:], paths))


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

    k = 35
    frequency_threshold = 3

    kmers = preprocess_pairs(input_reads, k, frequency_threshold)
    dbg = kmer_to_dbg(kmers)
    contigs = dbg_to_contigs(dbg)

    output_fn = args.output_file
    zip_fn = output_fn + '.zip'
    with open(output_fn, 'w') as output_file:
        output_file.write('>' + args.output_header + '\n')
        output_file.write('>ASSEMBLY\n')
        output_file.write('\n'.join(contigs))
    with zipfile.ZipFile(zip_fn, 'w') as myzip:
        myzip.write(output_fn)
