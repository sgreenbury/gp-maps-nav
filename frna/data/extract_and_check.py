#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import glob

df = pd.read_csv("frnadb_summary.zip")
lengths = np.arange(20, 141, 1)
df_sub = df[df["Length"].isin(lengths)].sort_values("Length")
(
    df_sub["Sequence"]
    .str.replace("T", "U")
    .to_frame()
    .to_csv("all_seqs_raw.csv", header=None, index=None)
)

df_from_src = pd.read_csv("all_seqs_raw.csv", header=None)
df_in_cs = pd.concat(
    [pd.read_csv(el, header=None) for el in glob.glob("cs_l*.txt")]
)

print(
    "Number not in intersection: {}".format(
        (~df_from_src[0].isin(df_in_cs[0])).sum()
        + (~df_in_cs[0].isin(df_from_src[0])).sum()
    )
)
print(
    "Number in intersection: {}".format(
        (df_from_src[0].isin(df_in_cs[0])).sum()
    )
)
