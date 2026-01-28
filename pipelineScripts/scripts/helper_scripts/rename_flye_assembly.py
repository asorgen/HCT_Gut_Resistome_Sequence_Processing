#!/usr/bin/env python3
from __future__ import print_function

import numpy as np # linear algebra
import pandas as pd # data processing, CSV file I/O (e.g. pd.read_csv)
import os
import sys

fasta = sys.argv[1]
info = sys.argv[2]

df = pd.read_csv(info, sep="\t")

for line in open(fasta):
    if line.startswith(">"): 
        line = line.strip()
        header = line[1:]
        # print(header)
        for index, row in df.iterrows():
            # print(f"Row {index}: A = {row['#seq_name']}, B = {row['length']}, C = {row['cov.']}")
            if row['#seq_name'] == header:
                print(f">{row['#seq_name']}_length_{row['length']}_cov_{row['cov.']}")
    else:
        print(line.strip())

    # for i in df["#seq_name"]:
    #     #     # print(i)
    #         if i == header: 
    #             print(i + " is the same as " + header)

        

