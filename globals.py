import datetime

import pandas as pd

print 'importing globals'

iter = 1

ratios, genome_seqs, anno, op_table, string_net, allSeqs_fname, all_bgFreqs, all_genes = \
    (None, None, None, None, None, None, None, None)

pynkey_code = {}

startTime = datetime.datetime.now()
print str(startTime)

clusters = {}

endTime = startTime

stats_df = pd.DataFrame()

