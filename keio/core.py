#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 12:26:50 2020

@author: Ravin Poudel
KEIO: A python module to process illumina reads for keio-collection type project.


"""

import os
import sys
import textdistance
import logging
import hashlib
import subprocess
import Bio
import gzip
import pickle
import nmslib
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger('keio.core')


def is_gzip(filename):
    try:
        with open(filename, "rb") as f:
            logging.info("check if %s is gzipped" % filename)
            return f.read(2) == b'\x1f\x8b'
    except IOError as e:
        logging.error("Could not open the file %s to determine if it was gzipped" % filename)
        raise e


def fq2fa(filelist, tempdir=None):
    """Saves a Fasta and from 1 or more fastq files (may be gzipped)

    Args:
        filelist (str): Genbank file to process

    Returns:
        None
    """
    try:
        fastpath = os.path.join(tempdir, "forward.fasta")
        with open(fastpath, "w") as f1:
            for file in filelist:
                if is_gzip(file):
                    with gzip.open(file, 'rt') as f:
                        records = SeqIO.parse(f, "fastq")
                        SeqIO.write(records, f1, "fasta")
                else:
                    with open(file, 'r') as f:
                        records = (SeqIO.parse(f, "fastq"))
                        SeqIO.write(records, f1, "fasta")
        return fastpath
    except Exception as e:
        print("An error occurred in input fastq file %s" % file)
        raise e

                    
def run_vsearch(mapping_fasta, reads_fasta, cluster_id=0.75, minseq_length=5, tempdir=None, threads=2):
    """ Returns mapping information
    Args:
            input1(str): barcodefile: fasta
            input2(str): reads: fastafile
            input3 (int): cluster_id
    Returns:
        output: file: vsearch output containing alignment positon and quality
    """
    try:
        out_info = 'query+target+ql+tl+id+tcov+qcov+ids+gaps+qrow+trow+id4+qilo+qihi+qstrand+tstrand'
        outputfile  = os.path.basename(mapping_fasta) + "__output.txt"
        fastpath = os.path.join(tempdir, outputfile)
        # print(reads_fasta)
        # print(mapping_fasta)
        # print(outputfile)
        # print(fastpath)
        parameters = ["vsearch", "--usearch_global", str(reads_fasta),
                      "--db", str(mapping_fasta), "--id", str(cluster_id),
                      "--minseqlength", str(minseq_length),
                      "--userfield", out_info,
                      "--strand", "both",
                      "--threads", str(threads),
                      "--userout", str(fastpath)]
        p0 = subprocess.run(parameters, stderr=subprocess.PIPE)
        print(p0.stderr.decode('utf-8'))
    except subprocess.CalledProcessError as e:
        print(str(e))


def parse_vsearch(file):
    """Parse vsearch file returning a dictionary of top hits for each primer and seqself.
    Args:
        input (str): file: the input file to be parsed
    Returns:
        (dict): A dictionary, e.g. {seq1: {primer1: {spos,epos,pmatch,tcov, gaps},...},...}
    """
    with open(file, 'r') as f:
        sdict = {}
        for line in f:
            ll = line.strip().split()
            qname = ll[0]
            tname = ll[1]
            pmatch = float(ll[4])
            tcov = float(ll[5])
            gaps = int(ll[8])
            spos = int(ll[12])
            epos = int(ll[13])
            qstrand = ll[14]
            tstrand = ll[15]
            td = {'spos': spos, 'epos': epos, 'pmatch': pmatch, 'tcov': tcov,
                  'gaps': gaps, 'qstrand': qstrand, 'tstrand': tstrand}
            if tcov >= 50:
                if gaps < 5:
                    if qname in sdict:
                        if tname in sdict[qname]:
                            if sdict[qname][tname]['pmatch'] < pmatch:
                                sdict[qname][tname] = td
                        else:
                            sdict[qname][tname] = td
                    else:
                        sdict[qname] = {}
                        sdict[qname][tname] = td
        return sdict




# filter the parsed vserach file based on the number of matching barcode type
def filter_vsearch(sdict, nhits):
    """Filter the parsed vserach file based on the number of matching barcode type"""
    outdict = {}
    for items in sdict.items():
        if len(items[1]) >= nhits:
            outdict[items[0]] = items[1]
    return outdict

# Now create a dictionary with start and end position for each barcode type


def get_randombarcode(keio_fasta, filter_vsearch_dict):
    """Create a dictionary with start and end position for each barcode type"""
    out_dict={}
    for seq_record in SeqIO.parse(keio_fasta, "fasta"):
        if seq_record.id in filter_vsearch_dict.keys():
            sequence = str(seq_record.seq)
            listkey = list(filter_vsearch_dict[seq_record.id].keys())
            a = filter_vsearch_dict[seq_record.id][listkey[0]]['spos'] # start position for Ffasta
            b = filter_vsearch_dict[seq_record.id][listkey[0]]['epos'] # start position for Ffasta
            c = filter_vsearch_dict[seq_record.id][listkey[1]]['spos']
            d = filter_vsearch_dict[seq_record.id][listkey[1]]['epos']
            up_constant = sorted(listkey)[1]
            down_constant = sorted(listkey)[0]
            ll = sorted([a, b, c, d])
            sp = ll[1]
            ep = ll[2]
            seq = sequence[sp:ep-1]
            if len(seq) >=18 and len(seq) <=22:
                out_dict[seq_record.id]= {"cutseq":seq, "up_constant":up_constant,"down_constant":down_constant}
    return out_dict
                

# Just get a fasta information for random barcode.
def randombarcode_fasta(get_randombarcode_dict):
    """ Retrive fasta from dictionary """
    barhash = []
    for k in get_randombarcode_dict.keys():
        rb = get_randombarcode_dict[k]['cutseq']
        record = SeqRecord(Seq(rb), id=k ,description="", name="")
        barhash.append(record)
    SeqIO.write(barhash, "rb.fasta", "fasta")


def cluster_db(rbfasta, threads=1, cluster_id=0.9, min_seqlength=10):
    """Runs Vsearch clustering to create a FASTA file of non-redundant sequences. Selects the most abundant sequence as the centroid
    Args:
        threads (int or str):the number of processor threads to use

    Returns:
            (file): uc file with cluster information
            (file): a centroid fasta file
    """
    try:
        centroid_fasta = rbfasta.split(".")[0] + "_centroid_representative_fasta"
        uc_file = rbfasta.split(".")[0] + ".uc"
        parameters0 = ["vsearch",
                      "--cluster_size", rbfasta,
                      "--id", str(cluster_id),
                      "--sizeout", "--sizeorder","--relabel",
                      "Cluster_",
                      "--centroids", centroid_fasta,
                      "--uc", uc_file,
                      "--strand", "both",
                      "--minseqlength", str(min_seqlength),
                      "--threads", str(threads)]
        p0 = subprocess.run(parameters0, stderr=subprocess.PIPE)
        print(p0.stderr.decode('utf-8'))
    except subprocess.CalledProcessError as e:
        print(str(e))
    except FileNotFoundError as f:
        print(str(f))


def mapR2clusterdb(fastafile, centroid_representative_fasta,cluster_id=0.9, min_seqlength=20):
    """Map reads and cluster centroid information"""
    try:
        uc_map = fastafile.split(".")[0] + ".cluster_table_mapping.uc"
        readfile = fastafile
        centroidfile = centroid_representative_fasta
        parameters1 = ["vsearch",
                      "--usearch_global", readfile,
                      "--db", centroidfile,
                      "--strand", "plus",
                      "--id", str(cluster_id),
                      "--uc", uc_map,
                      "--strand", "both",
                      "--minseqlength", str(min_seqlength)]
        p1 = subprocess.run(parameters1, stderr=subprocess.PIPE)
        print(p1.stderr.decode('utf-8'))
    except subprocess.CalledProcessError as e:
        print(str(e))
    except FileNotFoundError as f:
        print(str(f))

# use if NMSLIB
# index the reference list
def create_index(strings):
    """Create a nmslib index"""
    index = nmslib.init(space='leven',
                            dtype=nmslib.DistType.INT,
                            data_type=nmslib.DataType.OBJECT_AS_STRING,
                            method='small_world_rand')
    index.addDataPointBatch(strings)
    index.createIndex(print_progress=True)
    return index

# get knn in bactch mode for all query
def get_knns(index, vecs):
    return zip(*index.knnQueryBatch(vecs, k=1, num_threads=4)) # zip creates a tupple with first element as knn location and second element as distance

# Display the actual strings_KNN
def display_knn(ref_barcode_list, knn_ids_array):
    """Display the KNN neighbours"""
    #print("Query string:\n", query,"\n")
    print ("Pritning nearest neighbours:\n")
    for v in knn_ids_array.tolist():
        print(ref_barcode_list[v])


### filter based on distance
def filter_knn_dist(qdict, mindist):
    """Filter KNN based on distance"""
    outdict = {}
    for items in qdict.items():
         if items[1]['distance'] <= mindist:
            outdict[items[0]] = items[1]
    return outdict

def head_dict(dict, n):
    """Print first 'n ~ size' number of entries in a dictionary.

    Args:
        input (str): dict: the input dictionary to be parsed

        input (int): n: interger specifying the size of return dictionary

    Returns:
        print n entries for dictionary
    """
    k = 0
    for items in dict.items():
        if k < n:
            print("\n", items)
            k += 1
