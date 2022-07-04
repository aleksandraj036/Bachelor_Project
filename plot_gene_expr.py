from asyncio import subprocess
import os
import matplotlib.pyplot as plt
import pandas as pd
import plotly.express as px
import os
import numpy as np
import subprocess
import seaborn as sns
import plotly.io as pio

pio.kaleido.scope.mathjax = None #prevents from unnecessary error  popping up for pdf files

with open("OUT/gene_expression_all.csv", "w+") as genes:
        genes.write(f"Sample ID\tAll_genes\tExpressed_genes\n")
folder = "SALMON_OUT"

for root, dirs, files in os.walk(folder):
    for file_name in files:
        if "quant.sf" in file_name: 
            sample = root.split('/')[1]
            all_genes = sum(1 for line in open(f"{root}/{file_name}")) - 1
            dataframe = pd.read_csv(f"{root}/{file_name}",delimiter="\t")
            dataframe.to_csv(f"./OUT/salmon_{sample}.csv", encoding='utf-8', index = 0)
            data = pd.read_csv(f"./OUT/salmon_{sample}.csv" )
            df = pd.DataFrame(data)
            value = (df.iloc[:,3])
            no_zeros = value[(value > 0)]
            no_zeros = np.log10(no_zeros)
            df['zeros'] = no_zeros
            #creating density plot for every sample
            density = sns.displot(data = df, x= 'zeros',  kind="kde")
            density.set(xlabel = "Intervals", title=sample)
            density.figure.savefig(f"plots/density/{sample}.pdf", bbox_inches='tight') 
            expressed_genes = no_zeros.size
            with open("OUT/gene_expression_all.csv", "a") as genes:
                genes.write(f"{sample}\t{all_genes}\t{expressed_genes}\n")

#creating bar chart with gene expression %
with open("OUT/gene_expression_all.csv", "r") as g:
    dataframe2 = pd.read_csv(g, delimiter="\t")
    samples = dataframe2['Sample ID']
    all_g = dataframe2['All_genes']
    expressed = dataframe2['Expressed_genes']
    new = (expressed.astype(int)).divide(all_g.astype(int)) * 100
    dataframe2['Percentage'] = new
    un = px.bar(dataframe2, x="Sample ID", y="Percentage", title="Gene expression %")
    un.update_traces(marker_pattern_shape="\\", marker_color="#f48533")
    un.write_html("plots/density/pct_gene_expression.html")
    un.write_image('plots/density/pct_gene_expression.pdf')


