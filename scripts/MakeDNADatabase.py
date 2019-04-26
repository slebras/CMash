#! /usr/bin/env python
# This script will make a training database of hashes
import os
import sys
import math
import pickle
import requests
# The following is for ease of development (so I don't need to keep re-installing the tool)
try:
    from CMash import MinHash as MH
except ImportError:
    try:
        import MinHash as MH
    except ImportError:
        sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        from CMash import MinHash as MH
from multiprocessing import Pool  # Much faster without dummy (threading)
import multiprocessing
from itertools import *
import argparse
import khmer
import gzip
import shutil
from khmer.khmer_args import optimal_size\

# This function will make a single min hash sketch upon being given a file name, sketch size, prime, and k-mer size
def make_minhash(genome, max_h, prime, ksize):
    MHS = MH.CountEstimator(n=max_h, max_prime=prime, ksize=ksize, save_kmers='y', input_file_name=genome)
    # Just use HLL to estimate the number of kmers, no need to get exact count
    hll = khmer.HLLCounter(0.01, ksize)
    hll.consume_seqfile(genome)
    MHS._true_num_kmers = hll.estimate_cardinality()
    MHS.input_file_name = genome
    return MHS


# Unwrap for Python2.7
def make_minhash_star(arg):
    return make_minhash(*arg)


def stream_file(url):
    resp = requests.get(url, stream=True)
    if resp.status_code == requests.codes.ok:
        file_path = "fastas/"+ url.split("/")[-1]

        # stream file from response object
        with open(file_path, "ab") as fd:
            for chunk in resp.iter_content(chunk_size=2048):
                fd.write(chunk)

        # unzip the fasta file
        with gzip.open(file_path, "rb") as f_in:
            with open(file_path[:-6], "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)

        os.remove(file_path)
    else:
        raise Exception("Request for %s finished with error"%url)

    return file_path[:-6]


def unzip_file(file_path, temp_path):
    fasta_path = os.path.join(temp_path, file_path[:-6].split('/')[-1])
    with gzip.open(file_path, 'rb') as f_in:
        with open(fasta_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    return fasta_path

def main():
    parser = argparse.ArgumentParser(description="This script creates training/reference sketches for each FASTA/Q file"
                                    " listed in the input file.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-p', '--prime', help='Prime (for modding hashes)', default=9999999999971)
    parser.add_argument('-t', '--threads', type=int, help="Number of threads to use", default=multiprocessing.cpu_count())
    parser.add_argument('-n', '--num_hashes', type=int, help="Number of hashes to use.", default=500)
    parser.add_argument('-k', '--k_size', type=int, help="K-mer size", default=21)
    parser.add_argument('-i', '--intersect_nodegraph', action="store_true", \
                        help="Optional flag to export Nodegraph file (bloom filter) containing all k-mers in the" \
                             " training database. Saved in same location as out_file. This is to be used with QueryDNADatabase.py")
    # adding new parser argument to handle
    parser.add_argument('-t', '--temp_dir', type=str, help="temporary storage directory (define for continue flag)", default="./temp")
    parser.add_argument('-s', '--data_stream', action="store_true", \
                        help="Optional flag to define whether the input_files are urls to stream data instead of" \
                             " absolute paths to files.", default=False)
    parser.add_argument('-z', '--unzip_data', action="store_true", \
                        help="Optional flag to define whether the input_files are gzipped. if True, will unzip in " \
                             "chunks and delete unzipped fastas after use", default=False)
    parser.add_argument('-c', '--continue', action="store_true", \
                        help="Optional flag to define whether to continue sketching files defined in input file. " \
                             "Functionally, checks against the existing sketches in the temporary directory.", default=False)
    parser.add_argument('in_file', help="Input file: file containing (absolute) file names of training genomes.")
    parser.add_argument('out_file', help='Output training database/reference file (in HDF5 format)')
    args = parser.parse_args()
    num_threads = args.threads
    prime = args.prime  # taking hashes mod this prime
    ksize = args.k_size
    if ksize > 31:
        raise Exception("Unfortunately, ksize must be size 32 or smaller (due to khmer contraints). Please reduce the ksize or use MakeStreamingDNADatabase.py instead.")
    max_h = args.num_hashes
    input_file_names = os.path.abspath(args.in_file)
    if not os.path.exists(input_file_names):
        raise Exception("Input file %s does not exist." % input_file_names)
    out_file = os.path.abspath(args.out_file)
    if args.intersect_nodegraph is True:
        intersect_nodegraph_file = os.path.splitext(out_file)[0] + ".intersect.Nodegraph"
    else:
        intersect_nodegraph_file = None
    # create temporary dict if it doesn't exist
    if not os.path.isdir(args.temp_dir):
        os.mkdir(args.temp_dir)

    if args.unzip_data is True and args.data_stream is True:
        raise InputError("unzip_data and data_stream flags cannot both be specified.")

    if args.unzip_data is True or args.data_stream is True:
        chunk_size = 40
        with open(input_file_names, 'r') as fid:
            lines = fid.readlines()
        chunks = []
        for i in range(int(math.ceil(len(lines) / chunk_size))):
            if (i+1)*chunk_size > len(lines)-1:
                chunks[i*chunk_size:len(lines)]
            else:
                chunks[i*chunk_size:(i+1)*chunk_size]

    genome_sketches = []

    temp_path = args.temp_dir
    if args.unzip_data:
        if not os.path.isdir(os.path.join(temp_path, "fastas")):
            os.mkdir(os.path.join(temp_path, "fastas"))
        for idx, chunk in enumerate(chunks):
            print("Beginning download of chunk %i of %i"%(idx, len(chunks)))
            file_names = []
            for line in chunk:
                if not check_if_pickled(line):
                    f = unzip_file(line, os.path.join(temp_path, "fastas"))
                    file_names.append(f)

            if len(file_names) > 0:
                print("starting sketches")
                pool = Pool(processes=num_threads)
                curr_genome_sketches = pool.map(make_minhash_start, zip(file_names, repeat(max_h), repeat(prime), repeat(ksize)))
                genome_sketches += curr_genome_sketches

                pickle_sketches(curr_genome_sketches, temp_path)

                print("removing fasta files")
                for file_name in file_names:
                    os.remove(file_name)
            else:
                print("pickled files found, continuing...")

    # adding new
    elif args.data_stream is True:

        for idx, chunk in enumerate(chunks):
            print("Beginning download of chunk %i of %i"%(idx, len(chunks)))
            file_names = []
            for line in chunk:
                file = stream_file(line.strip())
                file_names.append(file)
            print("starting sketches")

            pool = Pool(processes=num_threads)
            curr_genome_sketches = pool.map(make_minhash_start, zip(file_names, repeat(max_h), repeat(prime), repeat(ksize)))
            genome_sketches += curr_genome_sketches

            print("removing fasta files")
            for file_name in file_names:
                os.remove(file_name)

    else:
        file_names = list()
        fid = open(input_file_names, 'r')
        for line in fid.readlines():
            line = line.strip()
            if not os.path.exists(line):
                raise Exception("Training genome %s does not exist." % line)
            file_names.append(line)
        fid.close()

        # Open the pool and make the sketches
        pool = Pool(processes=num_threads)
        genome_sketches = pool.map(make_minhash_star, zip(file_names, repeat(max_h), repeat(prime), repeat(ksize)))
    print("Beginning export to one HDF5 file")
    # Export all the sketches
    MH.export_multiple_to_single_hdf5(genome_sketches, out_file)

    # If requested, save all the k-mers into a big Nodegraph (unfortunately, need to pass through the data again since we
    # a-priori don't know how big of a table we need to make
    if intersect_nodegraph_file is not None:
        total_num_kmers = 0
        for sketch in genome_sketches:
            total_num_kmers += sketch._true_num_kmers
        res = optimal_size(total_num_kmers, fp_rate=0.001)
        intersect_nodegraph = khmer.Nodegraph(ksize, res.htable_size, res.num_htables)
        for file_name in file_names:
            intersect_nodegraph.consume_seqfile(file_name)
        intersect_nodegraph.save(intersect_nodegraph_file)

if __name__ == "__main__":
    main()
