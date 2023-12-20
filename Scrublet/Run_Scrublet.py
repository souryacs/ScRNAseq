# %matplotlib inline
import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
from optparse import OptionParser
import pandas as pd

parser = OptionParser()
parser.add_option("-m", 
                  "--matrix", 
                  dest="matrixfilename", 
                  help="Matrix filename")
parser.add_option("-g", 
                  "--gene", 
                  dest="genefilename", 
                  help="Gene filename")
parser.add_option("-o", 
                  "--out", 
                  dest="outfilename", 
                  help="Output filename")

(options, args) = parser.parse_args()

counts_matrix = scipy.io.mmread(options.matrixfilename).T.tocsc()
genes = np.array(scr.load_genes(options.genefilename, delimiter='\t', column=1))

print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
print('Number of genes in gene list: {}'.format(len(genes)))

scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets()

## dump the scrublet scores
df = pd.DataFrame({
    'doublet_score': scrub.doublet_scores_obs_,
    'predicted_doublet': scrub.predicted_doublets_
})
df.to_csv(options.outfilename, index=False)


