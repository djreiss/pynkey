import datetime

import pandas as pd

print 'importing globals'

startTime = datetime.datetime.now()
print str(startTime)

iter = 1

ratios, genome_seqs, anno, op_table, string_net, allSeqs_fname, all_bgFreqs, all_genes = \
    (None, None, None, None, None, None, None, None)

pynkey_code = {}

clusters = {}

endTime = startTime

stats_df = pd.DataFrame()


import logging

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s') #,
                    #filename='myapp.log',
                    #filemode='w')
#logging.debug('A debug message')
#logging.info('Some information')
#logging.warning('A shot across the bows')
