import datetime

import pandas as pd
import numpy as np

import init
import params
##from Bicluster import bicluster

print 'importing globals'

iter = 1

ratios, genome_seqs, anno, op_table, string_net, allSeqs_fname, all_bgFreqs, all_genes = \
    init.pynkey_init(params.organism, params.k_clust, params.ratios_file)

# Save all pynkey code for safe keeping
pynkey_code = {}
try:
    pynkey_code = init.load_pynkey_code()
except:
    print 'Must have "import"ed this module rather than run from main.py.'
    
startTime = datetime.datetime.now()
print str(startTime)

clusters = init.init_biclusters( ratios, params.k_clust, 'kmeans+random' );
        
##counts_g = bicluster.get_all_cluster_row_counts( clusters, all_genes )

print 'DONE WITH INITIALIZATION!'
endTime = datetime.datetime.now()
print str(endTime)
print str(endTime - startTime) + ' seconds since initialization'
        
stats_df = pd.DataFrame()

