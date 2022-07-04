import subprocess
import os
import sys
import time

'''
The root folder must contain:
1. File with rRNA for analyzed organism
2. Folder "samples" with all the samples used for analysis
3. File with genome of analyzed organism in FASTA format
4. File with known splice sites in GTF format.
5. File with coding transcripts
6. File with non-coding transcripts
7. File differential_expression.R
8  Scripts: plot mapping.py, plot gene expr.py, plot qc.py, prepare for deseq salmon.py, differential expression.R
9. File expression_samples.txt <- needs to be manually created by the user (description in thesis paper)
10. File with transcript<->gene annotation can be downloaded from ENSEMBL using BioMart
'''
# run as 'python3 main.py /home/username/rnaseq/samples'
#                script_name   current directory with samples


# Check if all necessary programs are installed*********************************************

programs = ['hisat', 'stringtie', 'trimmomatic']
for p in programs:
    if not os.path.isdir(f'{os.getcwd()}/{p}'):
        print(f'Installing {p}')
        subprocess.run(f'wget http://rhesus.amu.edu.pl/share/programs/{p}.zip', shell=True)
        subprocess.run(f'unzip {p}.zip', shell=True)

#installing Salmon
subprocess.run(f'wget https://github.com/COMBINE-lab/salmon/releases/download/v0.9.1/Salmon-0.9.1_linux_x86_64.tar.gz',shell=True)
subprocess.run(f'tar -xvf Salmon-0.9.1_linux_x86_64.tar.gz',shell=True)

print('All necessary programs are installed')

# ******************************************************************************************

# Take folder name from the user and iterate over samples in this folder********************
folder = sys.argv[1]  # user needs to provide name of the folder where are all the samples why calling the file


# FILL OUT
# ******************************************************************************************
species = 'human' 

sample_naming = 1 #if sample naming is _1 and _2 should be equal to 1, if it is _R1 and _R2 should be equal to 0

genome_fasta = 'Homo_sapiens.GRCh38.dna.primary_assembly.fa'

genome_gtf = 'Homo_sapiens.GRCh38.106.gtf'

genome = 'Homo_sapiens.GRCh38'

folder_in_index = 'index\chr22'

group_1 = 'control'

group_2 = 'case'
# ******************************************************************************************


sample_1 = ''
sample_2 = ''
sample_id = ''


def qc(sample_id):
    if sample_naming:
        sample_1 = sample_id + "_1.fastq"
        sample_2 = sample_id + "_2.fastq"
    else:
        sample_1 = sample_id + "_R1.fastq"
        sample_2 = sample_id + "_R2.fastq"

    quality_control = [
            f'java -Xms4g -Xmx4g -jar trimmomatic/trimmomatic.jar PE -threads 1 -phred33 samples/{sample_1} samples/{sample_2} '
            f'TRIMMED/{sample_id}_trimmomatic_R1_trimmed.fastq /dev/null TRIMMED/{sample_id}_trimmomatic_R2_trimmed.fastq '
            f'/dev/null ILLUMINACLIP:trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:50 2> '
            f'STATUS/{sample_id}_trimmomatic.log',
            f'(bowtie2 -t -p 4 -X 1000 -1 TRIMMED/{sample_id}_trimmomatic_R1_trimmed.fastq -2 '
            f'TRIMMED/{sample_id}_trimmomatic_R2_trimmed.fastq -x index/{species}_rRNA --fast --un-conc {sample_id}.fastq) >/dev/null '
            f'2>bowtie2/{sample_id}_bowtie.log',

            f'mv samples/{sample_1} TRIMMED/{sample_id}_clean_R1.fastq',

            f'mv samples/{sample_2} TRIMMED/{sample_id}_clean_R2.fastq'
        ]
    return quality_control

def genome_assembly(): 
    genome_assembly = [
            f'chmod u+x -R hisat',
            f'hisat/hisat-build {genome_fasta} {folder_in_index}', 
            f'python hisat/extract_splice_sites.py {genome_gtf} > {species}_splice_sites.txt',
    ]
    return genome_assembly

def mapping(sample_id):
    mapping = [
            f'chmod u+x -R hisat',
            f'hisat/hisat -q -p 4 -X 1000 --time --phred33 --rna-strandness RF --known-splicesite-infile '
            f'{species}_splice_sites.txt --novel-splicesite-outfile novel.splice_sites.txt -x {folder_in_index} -1 '
            f'TRIMMED/{sample_id}_clean_R1.fastq -2 TRIMMED/{sample_id}_clean_R2.fastq -S {sample_id}.sam 2> '
            f'STATUS/{sample_id}_hisat.txt',
            f'samtools view -bS {sample_id}.sam > {sample_id}.bam',

            f'rm {sample_id}.sam',

            f'samtools sort {sample_id}.bam -o {sample_id}.sorted.bam',

            f'samtools index {sample_id}.sorted.bam',

            f'stringtie {sample_id}.sorted.bam -o {sample_id}.gtf -p 1 -G {genome_gtf} -A abundance.txt',
            f'cuffcompare -r {genome_gtf} -R -o {sample_id} -C -G {sample_id}.gtf',
            f'python counts_stringtie.py {sample_id}.combined.gtf' 
        ]
    return mapping


def salmon_index():
    create_index = [
            f'cat {genome}.cdna.all.fa {genome}.ncrna.fa > {species}_transcriptome.fasta',
            f'Salmon-latest_linux_x86_64/bin/salmon index -p 2 -t {species}_transcriptome.fasta -i index/{species}_transcriptome --type quasi -k 31' #tworzenie indeksu
    ]
    return create_index

def gene_expr(sample_id):
    expression_estimation = [
            f'Salmon-latest_linux_x86_64/bin/salmon quant -p 4 --useVBOpt -i index/{species}_transcriptome -l IU -1 TRIMMED/{sample_id}_clean_R1.fastq -2 TRIMMED/{sample_id}_clean_R2.fastq -o SALMON_OUT/{sample_id}'
    ]
    return expression_estimation


paired_end = {} #dictionary that will contain ID of all the samples


for root, dirs, files in os.walk(folder):
        for file_name in files:
            if "_" in file_name: #paired end
                paired = file_name.split('_')[0]
                if paired not in paired_end:
                    paired_end[paired] = 1


folders = ['OUT', 'TRIMMED', 'STATUS', 'index', 'bowtie2', 'SALMON_OUT', 'plots', 'plots/differential', 'plots/density']

#create necessary folders
for folder_name in folders:
    subprocess.run(f'mkdir {folder_name}', shell=True)

#build bowtie index
subprocess.run(f'bowtie2-build {species}_rRNA.fasta index/{species}_rRNA',  shell=True)

#run RNA-Seq
for step in [qc, genome_assembly, mapping, salmon_index, gene_expr]:
    #st = time.time()   #ALL COMMENTED LINES IN THIS PART WERE USED FOR MEASURING EXECUTION TIME
    if (step.__name__ == "genome_assembly") or (step.__name__ == "salmon_index"):
        for command in step(): 
            subprocess.run(command, shell=True)
        #et = time.time()
        #with open('time.txt', "a+" ) as t:
        #    t.write(step.__name__ + "   ")
        #    t.write(str(et - st) + "\n")
        continue
    for sample in paired_end:
        for command in step(sample):
            subprocess.run(command, shell=True)
    subprocess.run(f'python3 plot_{step.__name__}.py', shell=True)
    #et = time.time()
    #with open('time.txt', "a+" ) as t:
    #    t.write(step.__name__ + "   ")
    #    t.write(str(et - st) + "\n")
    


#preparing input for DESeq2
#st = time.time() #AGAIN USED TO MEASURE TIME
subprocess.run(f'python3 prepare_for_deseq_salmon.py expression_samples.txt {group_1} {group_2} SALMON_OUT ensembl_data_{species}.txt {group_2}_vs_{group_1}.txt',shell=True)


#differential expression using DESeq2
os.system(f'Rscript differential_expression.R {group_2}_vs_{group_1}.txt')

#et = time.time() #TIME
#with open('time.txt', "a+" ) as t:
#    t.write("differential_expression" + "   ")
#    t.write(str(et - st) + "\n")

