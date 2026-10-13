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
    print('bnumber\tdbname\tdbid')
    for f in record.features:
        if f.type=='CDS':
            #print(f)
            f_loc = f.location
            f_xref = f.qualifiers['db_xref']
            f_id = f.qualifiers['locus_tag']
            f_desc = f.qualifiers['product']
            for i in f_id:
                for j in f_xref:
                    (dbsource, dbid) = j.split(':')
                    dbsource = 'UniProtKB' if dbsource.startswith('UniProtKB') else dbsource
                    print(f'{i}\t{dbsource}\t{dbid}')