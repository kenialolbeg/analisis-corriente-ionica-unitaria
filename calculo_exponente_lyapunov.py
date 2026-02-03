###Sección para calcular Lyap de K_channel ECKMANN
import numpy as np
import nolds
import pandas as pd
colnames=['K']
df = pd.read_csv("K_channel_normalized.csv",names=colnames, header=None, sep=',')

lyap = nolds.lyap_e(df.K[1:100000], emb_dim=5, matrix_dim=5, min_nb=None,min_tsep=9000, tau=1, debug_plot=False)
print(lyap)