import sys

log_file=open(sys.argv[1]).read().split('\n')[:-1]
num = 1
print("gene_callers_id\tcontig\tstart\tstop\tdirection\tpartial\tsource\tversion")
for gene in log_file[:10]:
	gene_name = gene.split('\t')[0]
	gene_pos = gene.split('\t')[1].split('_')
	start = gene_pos[-3]
	stop = gene_pos[-2]
	if gene_pos[-1] == '+':
		direction = 'f'
	else:
		direction = 'r'
	contig = '_'.join(gene_pos[:-3])
	print(str(num)+ "\t" + contig + "\t" + start + "\t" + stop + "\t" + direction + "\t0\tFragGeneScan\t1.3")
	num += 1
