#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 12:31:42 2020

@author: Ravin Poudel
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
import gzip
import keio

#######


def test_run_vsearch():
	keio.run_vsearch(mapping_fasta ="test/test_data/upstream.fasta", tempdir = "test/test_data/")
	