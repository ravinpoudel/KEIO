#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 12:31:42 2020

@author: admin
"""
import os
import sys
import logging
import argparse
import tempfile
import shutil
import textdistance
import hashlib
import subprocess
import Bio
import pickle
import nmslib
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
from Bio.SeqRecord import SeqRecord



import keio


##subprocess 

# reformat.sh in=reads.fastq out=reads.fasta

############


def myparser():
    parser = argparse.ArgumentParser(description='Mapping inline barcodes to the fasta file')
    parser.add_argument('--fastq', '-f', nargs='+', type=str, required=True, help='input fastq file')
    parser.add_argument('--upstreamFasta', '-uf', type=str, required=True, help='A upstreamFasta file')
    parser.add_argument('--downstreamrcFasta', '-drcf', type=str, required=True, help='A downstreamFasta file')
    parser.add_argument('--threads', help='The number of cpu threads to use', type=int, default=2)
    #parser.add_argument('--log', help="Log file", default="keio.log")
    parser.add_argument('--tempdir', help='The temp file directory', default=None)
    parser.add_argument('--keeptemp' ,help="Should intermediate files be kept?", action='store_true')
    return parser



def _logger_setup(logfile):
    """Set up logging to a logfile and the terminal standard out.

    Args:
        logfile (str): Log file

    """
    try:
        logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                            datefmt='%m-%d %H:%M',
                            filename=logfile,
                            filemode='w')
        # define a Handler which writes INFO messages or higher to the sys.stderr
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)
        # set a format which is simpler for console use
        formatter = logging.Formatter('%(asctime)s: %(levelname)-8s %(message)s')
        # tell the handler to use this format
        console.setFormatter(formatter)
        # add the handler to the root logger
        logging.getLogger('').addHandler(console)
    except Exception as e:
        print("An error occurred setting up logging")
        raise e
        
def main(args=None):
    """Run The complete Keio workflow.

    """
    # Set up logging
    parser = myparser()
    if not args:
        args = parser.parse_args()

    basename = os.path.basename(str(args.fastq)) 
    logfilename = basename.split(".")[0] + ".log"
    _logger_setup(logfilename)

    try:

        if args.tempdir:
            if not os.path.exists(args.tempdir):
                logging.warning("Specified location for tempfile (%s) does not \
                                 exist, using default location." % tempdir )
                tempdir = tempfile.mkdtemp(prefix='keio_')
        else:
            tempdir = tempfile.mkdtemp(prefix='keio_', dir=args.tempdir)
            
        logging.info("Temp directory is: %s" % (tempdir))
        logging.info("Writing fasta file from genbank file(s)")
        fastapath = keio.fq2fa(args.fastq, tempdir=tempdir)
        print(fastapath)

        # upstream
        logging.info("Mapping upstream barcode -- Number of threads used: %d", (args.threads))
        keio.core.run_vsearch(mapping_fasta=args.upstreamFasta, reads_fasta = fastapath, cluster_id=0.95, minseq_length=8,  threads=args.threads, tempdir=tempdir)
        
        # downstream
        logging.info("Mapping downstream barcode -- Number of threads used: %d", (args.threads))
        keio.run_vsearch(mapping_fasta=args.downstreamrcFasta, reads_fasta = fastapath, cluster_id=0.90, minseq_length=8, threads=args.threads, tempdir=tempdir)

        # NEed to move this processing to core.
        f_list = [x for x in os.listdir(tempdir) if x.endswith(".txt")]
        
        outfilename = tempdir + '/PlateMappingReads_combinedfile.txt'

        logging.info("Combinging upstream and downstrem output at: %s", (outfilename))

        with open(outfilename, 'w') as outfile:
            for fname in f_list:
                file = os.path.join(tempdir, fname)
                with open(file) as infile:
                    for line in infile:
                        outfile.write(line)
        ###################
        
        logging.info("Parsing vsearch records")
        datplus = keio.parse_vsearch(outfilename)
        logging.info("Number of records: %d", (len(datplus)))
        
        logging.info("Filtering records")
        datplus_filter = keio.filter_vsearch(datplus, nhits=2)
        logging.info("Number of records: %d", (len(datplus_filter)))
        #print(len(datplus_filter))
        #keio.head_dict(datplus_filter, 5)


        basename = os.path.basename(str(args.fastq)) 
        outputfile = basename.split(".")[0] + ".csv"
        
        PlateMapping_rbdict = keio.get_randombarcode(fastapath, datplus_filter)
        #keio.head_dict(PlateMapping_rbdict, 5)
        
        df_PlateMapping_rbdict = pd.DataFrame.from_dict(PlateMapping_rbdict ,orient='index')
        df_PlateMapping_rbdict_ri = df_PlateMapping_rbdict.rename_axis('SeqID').reset_index()

        logging.info("Writing PlateMappingReads as csv file: %s", outputfile)
        df_PlateMapping_rbdict_ri.to_csv(outputfile, index=False)

    except Exception as e:
        logging.error("Keio terminated with errors. See the log file for details.")
        logging.error(e)
        raise SystemExit(1)
    finally:
        try:
            if not args.keeptemp:
                shutil.rmtree(tempdir)
        except UnboundLocalError:
            raise SystemExit(1)
        except AttributeError:
            raise SystemExit(1)

 

