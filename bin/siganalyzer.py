#!/usr/bin/env python

import sys
import pandas as pd
import numpy as np

sample=sys.argv[1]
path_script=sys.argv[2]
output_dir=sys.argv[3]

path_sig = output_dir + "nmf_output.h5"
path_cosmic = path_script + "/cosmic_ref/COSMIC_v2_SBS_GRCh38.txt"

H = pd.read_hdf(path_sig, 'H')
H = H.iloc[: , :-3]

matrix_samples = pd.read_csv(sample, sep='\t', index_col=0)
cosmic2 = pd.read_csv(path_cosmic, sep='\t', index_col=0)

df_columns = list(cosmic2.columns)
df_index = list(matrix_samples.columns)

zero_df = pd.DataFrame(0, index = df_index, columns = df_columns)
print(zero_df)

H_columns_formatted = [i.split('-')[1] for i in list(H.columns)]
H_columns_formatted = [i.replace(" ", "_") for i in H_columns_formatted]

H.columns = H_columns_formatted 

print(H.sort_index(axis=1))


#cols = df.columns.union(df2.columns)
