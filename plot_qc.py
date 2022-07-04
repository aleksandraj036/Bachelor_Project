import os
import pandas as pd
import plotly.express as px
import subprocess
import plotly.io as pio



with open(os.path.abspath(f"./OUT/samples.txt"), 'w+') as out:
    out.write("Sample ID" + "\t" +  "Raw" + "\t" "After QC" + "\t" +   "After rRNA")
    out.write('\n')


for root, dirs, files in os.walk('STATUS'):
    for file_name in files:
        if 'log' in file_name:
            current_file = os.path.join(root, file_name)
            with open(os.path.abspath(f"./{current_file}")) as log:
                num = [int(s) for s in ((log.readlines()[4]).lstrip(' ')).split() if s.isdigit()][:2]
                sample = file_name.split('_')[0]
                with open(os.path.abspath(f"./OUT/samples.txt"), 'a') as out:
                    out.write(sample + '\t')
                    for n in num:
                        out.write(str(n) + '\t')
                    with open(os.path.abspath(f"./bowtie2/{sample}_bowtie.log"), 'r') as rrna:
                        after_rrna = [int(s) for s in ((rrna.readlines()[6]).lstrip(' ')).split() if s.isdigit()][:2]
                        out.write(str(after_rrna[0]))
                        out.write('\n')



dataframe = pd.read_csv(f"./OUT/samples.txt",delimiter="\t")

dataframe.to_csv(f"./OUT/samples.tsv", encoding='utf-8', index = 0)

data = pd.read_csv(f"./OUT/samples.tsv" )

df = pd.DataFrame(data)

fig = px.bar(df, x="Sample ID", y=["Raw",  "After QC",  "After rRNA"], title="Quality Control Statistics").update_layout(yaxis_title="Number of reads")


pio.kaleido.scope.mathjax = None

#creating stacked bar chart
fig.write_html("plots/quality/samples_qc.html")
fig.write_image('plots/quality/samples_qc.pdf')
