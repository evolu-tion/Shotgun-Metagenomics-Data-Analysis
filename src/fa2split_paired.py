import sys, os, stat
import gzip
from optparse import OptionParser

# Parse Command Line
usage = """
Description: tool for split single paired-end file to fasta paired files
"""

parser = OptionParser(usage=usage)
parser.add_option("-i", "--input_single_fasta", dest="single_fasta",
	default=None, metavar="FILE",
	help="single paired-end fasta file")
parser.add_option("-o", "--output", dest="dir_out",
	default=None, metavar="FILE",
	help="Output fasta file name for example \'XXX\' for XXX_R1.paired.fa, XXX_R2.paired.fa, XXX.unpaired.fa")

options,args = parser.parse_args()

if not options.single_fasta:
	sys.exit("Missing input signle fasta file, -i xxx.fasta")
if not options.dir_out:
	sys.exit("Missing output directory, -o xxx")

dir_single_fasta = os.path.abspath(options.single_fasta)
dir_out_fa = options.dir_out

w_file_R1 = open(dir_out_fa+"_R1.paired.fa", 'w')
w_file_R2 = open(dir_out_fa+"_R2.paired.fa", 'w')

single_fasta_f = open(dir_single_fasta, 'r').read().split('\n')
single_fasta_f = single_fasta_f[:-1]

print("Begin to convert single fasta of "+ os.path.split(dir_out_fa)[1] +" with "+ str(len(single_fasta_f)/4) + " sequences.")

for i in range(0,len(single_fasta_f),4):
	if (single_fasta_f[i][:-2] == single_fasta_f[i+2][:-2]):
		w_file_R1.write(single_fasta_f[i][:-2] + "\n" + single_fasta_f[i+1] + "\n")
		w_file_R2.write(single_fasta_f[i+2][:-2] + "\n" + single_fasta_f[i+3] + "\n")
	else:
		sys.exit("Error of signle fasta file")

