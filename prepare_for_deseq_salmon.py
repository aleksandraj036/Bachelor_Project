import sys

samples_file = sys.argv[1]
control_features = sys.argv[2]
treatment_features = sys.argv[3]
salmon_dir = sys.argv[4]
ensembl_data = sys.argv[5]
output = sys.argv[6]


f = open(samples_file)
out = open(output, 'w')

control_samples, treatment_samples = [], []

with open(samples_file) as f:
	for line in f.readlines():
		if control_features in line:
			control_samples.append(line.split(None, 2)[0])
		else:
			treatment_samples.append(line.split(None, 2)[0])


print("Control:", control_features)
print("Treatment:", treatment_features)
print("Control samples: ", control_samples)
print("Treatment samples: ", treatment_samples)
	

	
dict_control, dict_treatment = {}, {}
for sample in control_samples:
	try:
		f.close()
	except:
		pass
	f = open(salmon_dir + '/' + sample + '/quant.sf')
	f.readline()
	for line in f:
		line = line.strip().split('\t')
		transcript = line[0].split('.')[0]
		expected_count = float(line[-1])
		if transcript in dict_control:
			dict_control[transcript].append(expected_count)
		else:
			dict_control[transcript] = [expected_count]


for sample in treatment_samples:
	try:
		f.close()
	except:
		pass
	f = open(salmon_dir + '/' + sample + '/quant.sf')
	f.readline()
	for line in f:
		line = line.strip().split('\t')
		transcript = line[0].split('.')[0]
		expected_count = float(line[-1])
		if transcript in dict_treatment:
			dict_treatment[transcript].append(expected_count)
		else:
			dict_treatment[transcript] = [expected_count]
			
	 

# getting gene-level expressions

dict_gene_transcript = {}
dict_transcript_gene = {}
for line in open(ensembl_data):
	line = line.split('\t')
	gene = line[0]
	transcript = line[1]
	if gene in dict_gene_transcript:
		dict_gene_transcript[gene].append(transcript)
	else:
		dict_gene_transcript[gene] = [transcript]
	dict_transcript_gene[transcript] = gene
	

dict_control_gene, dict_treatment_gene = {}, {}
for transcript in dict_control.keys():
	gene = dict_transcript_gene[transcript]
	expressions_transcript = dict_control[transcript]
	if gene in dict_control_gene:
		expressions_gene = dict_control_gene[gene]
		expressions_new = []
		for i in range(len(expressions_transcript)):
			expressions_new.append(expressions_transcript[i] + expressions_gene[i])
		dict_control_gene[gene] = expressions_new
	else:
		dict_control_gene[gene] = expressions_transcript
	
for transcript in dict_treatment.keys():
	gene = dict_transcript_gene[transcript]
	expressions_transcript = dict_treatment[transcript]
	if gene in dict_treatment_gene:
		expressions_gene = dict_treatment_gene[gene]
		expressions_new = []
		for i in range(len(expressions_transcript)):
			expressions_new.append(expressions_transcript[i] + expressions_gene[i])
		dict_treatment_gene[gene] = expressions_new
	else:
		dict_treatment_gene[gene] = expressions_transcript
	
	 
write = 'gene'
for i in range(len(control_samples)):
	write += '\tcontrol'  
for i in range(len(treatment_samples)):
	write += '\ttest'
	
out.write(write + '\n')     
for gene in dict_control_gene.keys():
	control_expressions = dict_control_gene[gene]
	treatment_expressions = dict_treatment_gene[gene]
	write = gene
	for control_expression in control_expressions:
		write += '\t' + str(int(round(control_expression)))
	for treatment_expression in treatment_expressions:
		write += '\t' + str(int(round(treatment_expression)))   
	write += '\n'
	out.write(write)
	