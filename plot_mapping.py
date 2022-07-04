from turtle import color
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import os
import re
import subprocess
import plotly.io as pio


with open(os.path.abspath(f"./OUT/hisat.txt"), 'w+') as out:
    out.write("Sample ID" + "\t" +  "Unmapped" + "\t" "Uniquely mapped" + "\t" +   "Multi-mapped" + "\t" + "%Unmapped" + "\t" "%Uniquely mapped" + "\t" +   "%Multi-mapped")
    out.write('\n')


for root, dirs, files in os.walk('STATUS'):
    for file_name in files:
        if 'hisat' in file_name:
            current_file = os.path.join(root, file_name)
            with open(os.path.abspath(f"./{current_file}")) as hisat:
                lines = hisat.readlines()
                unmapped_lst = re.findall(r'\d+(?:\.\d+)?', lines[5])
                unique_lst = re.findall(r'\d+(?:\.\d+)?', lines[6])
                multi_lst = re.findall(r'\d+(?:\.\d+)?', lines[7])
                sample = file_name.split('_')[0]
                with open(os.path.abspath(f"./OUT/hisat.txt"), 'a') as out:
                    out.write(sample + '\t' + unmapped_lst[0] +'\t' + unique_lst[0] +'\t' + multi_lst[0] +'\t' + unmapped_lst[1] +'\t' + unique_lst[1] +'\t' + multi_lst[1] + '\n')



dataframe = pd.read_csv(f"./OUT/hisat.txt",delimiter="\t")
dataframe.to_csv(f"./OUT/hisat.tsv", encoding='utf-8', index = 0)
data = pd.read_csv(f"./OUT/hisat.tsv" )
df = pd.DataFrame(data)



id = list(df.iloc[:, 0])

unmapped = list(df.iloc[:, 1])

uniq = list(df.iloc[:, 2])

multi = list(df.iloc[:, 3]) 
p_unmapped = list(df.iloc[:, 4])

p_uniq = list(df.iloc[:, 5])

p_multi = list(df.iloc[:, 6])


pio.kaleido.scope.mathjax = None
#unmapped
un = px.bar(df, x="Sample ID", y="Unmapped", title="Unmapped")
un.update_traces(marker_color="#1d73aa")
un.write_html("plots/mapping/unmapped.html")
un.write_image('plots/mapping/unmapped.pdf')

#uniquely
uni = px.bar(df, x="Sample ID", y= "Uniquely mapped" , title="Uniquely mapped")
uni.update_traces(marker_color="#c44441")
uni.write_html("plots/mapping/uniquely_mapped.html")
uni.write_image('plots/mapping/uniquely_mapped.pdf')
#multimapped
mlt = px.bar(df, x="Sample ID", y="Multi-mapped" , title="Multi-mapped")
mlt.update_traces(marker_color="#8cbc4f")
mlt.write_html("plots/mapping/multi-mapped.html")
mlt.write_image('plots/mapping/multi-mapped.pdf')

#unmapped %
pun = px.bar(df, x="Sample ID", y="%Unmapped", title="Percentage of unmapped")
pun.update_traces(marker_pattern_shape="\\", marker_color="#1d73aa")
pun.write_html("plots/mapping/pct_unmapped.html")
pun.write_image('plots/mapping/pct_unmapped.pdf')
#uniquely %
puni = px.bar(df, x="Sample ID", y= "%Uniquely mapped" , title="Percentage of uniquely mapped")
puni.update_traces(marker_pattern_shape="\\", marker_color="#c44441")
puni.write_html("plots/mapping/pct_uniquely_mapped.html")
puni.write_image('plots/mapping/pct_uniquely_mapped.pdf')
#multimapped %
pmlt = px.bar(df, x="Sample ID", y="%Multi-mapped" , title="Percentage of multi-mapped")
pmlt.update_traces(marker_pattern_shape="\\", marker_color="#8cbc4f")
pmlt.write_html("plots/mapping/pct_multi-mapped.html")
pmlt.write_image('plots/mapping/pct_multi-mapped.pdf')