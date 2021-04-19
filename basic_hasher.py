import sys
import argparse
import numpy as np
import time
import zipfile
import textwrap
from collections import defaultdict, Counter, OrderedDict


def parse_reads_file(reads_fn):
    """
    :param reads_fn: the file containing all of the reads
    :return: outputs a list of all paired-end reads

    HINT: This might not work well if the number of reads is too large to handle in memory
    """
    try:
        with open(reads_fn, 'r') as rFile:
            print("Parsing Reads")
            first_line = True
            count = 0
            all_reads = []
            for line in rFile:
                count += 1
                # if count % 1000 == 0:
                #     print(count, " reads done")
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


def create_ref_dict(reference, k):
    d = OrderedDict()
    for i in range(len(reference)-k):
        seq = reference[i:i+k]
        d.setdefault(seq, []).append(i)
    return d


def snp(reads, reference, k):
    genome = create_ref_dict(reference, k)

    def check_snp(read, genome_seq):
        index = None
        mismatches = 0
        for i in range(len(read)):
            if read[i] != genome_seq[i]:
                mismatches += 1
                index = i
        if not index:
            return None
        if 3 > mismatches:
            return[genome_seq[index], read[index], index]
        return None

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
            if not index_array:
                continue
            for i in index_array:
                start = i - first_match_index
                genome_seq = reference[start: start+read_length]
                x = check_snp(read, genome_seq)
                if x is not None:
                    snps.append([x[0], x[1], start+x[2]])
    snps.sort(key=lambda x: x[2])
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
    return snps_new


def detect_ins(read, ref, ref_start):  # read is longer than ref
    result = []
    k = len(read)
    ref_len = len(ref)
    for i in range(k-ref_len+1):
        if ref == read[i:i+ref_len]:
            if i > 0:
                result.append((read[:i], ref_start))
            elif i + ref_len < k:
                result.append((read[i+ref_len:], ref_start+i+ref_len))
            break
    return result


def detect_del(read, ref, ref_start):  # read is shorter than ref
    result = []
    k = len(read)
    ref_len = len(ref)
    for i in range(ref_len-k+1):
        if read == ref[i:i+k]:
            if i > 0:
                result.append((ref[:i], ref_start))
            elif i + k < ref_len:
                result.append((ref[i+k:], ref_start+i+k))
            break
    return result


def indel(reads, reference, k):
    genome = create_ref_dict(reference, k)
    ins = []
    dels = []
    limit = 3
    for read in reads:
        first = read[:k]
        read_mid = read[k:2*k]
        last = read[2*k:3*k]

        first_pos = genome.get(first)
        last_pos = genome.get(last)
        if first_pos and last_pos:
            for x_0 in first_pos:
                x = x_0 + k
                for y in last_pos:
                    if x+k-limit < y and y < x+k:  # insert
                        ins_strings = detect_ins(
                            read_mid, reference[x:y], x)
                        ins.extend(ins_strings)
                    elif x+k < y and y < x+k+limit:  # delete
                        del_strings = detect_del(
                            read_mid, reference[x:y], x)
                        dels.extend(del_strings)

    ins_final = list(set(ins))
    ins_final.sort(key=lambda x: x[1])
    dels_final = list(set(dels))
    dels_final.sort(key=lambda x: x[1])
    return (ins_final, dels_final)


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
    reads = []
    for pair in input_reads:
        reads.append(pair[0])
        reads.append(pair[1])

    k_snp = 25
    k_indel = 16

    snps = snp(reads, reference, k_snp)
    indels = indel(reads, reference, k_indel)
    insertions = indels[0]
    deletions = indels[1]

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
