#!/bin/env python
# conda install python3 biopython requests
# ./gbff.to.bnumber.GeneID.py > mapping.bnumber_ncbi.tsv

import gzip
import re, sys
import pycparser
from Bio import SeqIO, Entrez


gbff = sys.argv[1]# './downloads/GCF_000005845.2.gbff.gz'

with gzip.open(gbff, "rt") as handle:
    record = SeqIO.read(handle, "genbank")
    rank = 0
    print('bnumber\trank\tstrand\tbegin\tend')
    for f in record.features:
        if f.type=='CDS':
            rank +=1
            f_loc = str(f.location)
            m = re.search("\[(\d+):(\d+)\]\((.)\)",f_loc)
            f_xref = f.qualifiers['db_xref']
            f_id = f.qualifiers['locus_tag']
            f_desc = f.qualifiers['product']
            for i in f_id:
                print(f'{i}\t{rank}\t{m[3]}\t{m[1]}\t{m[2]}')
               
