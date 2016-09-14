import sys, os, stat
import gzip
from optparse import OptionParser

# Parse Command Line
usage = """
Description: tool for merage fastq file and convert to single fasta file
"""

parser = OptionParser(usage=usage)
parser.add_option("-1", "--fq1", dest="dir_fastq_1",
	default=None, metavar="FILE",
	help="Forward paired-end fastq file")
parser.add_option("-2", "--fq2", dest="dir_fastq_2",
	default=None, metavar="FILE",
	help="Backward paired-end fastq file")
parser.add_option("-o", "--output", dest="dir_out",
	default=None, metavar="FILE",
	help="Output fasta file name")

options,args = parser.parse_args()

if not options.dir_fastq_1:
	sys.exit("Missing input fastq 1 file, -1 xxx.fastq")
if not os.path.isfile(options.dir_fastq_1):
	sys.exit("Missing input fastq 1 file, %r" % options.dir_fastq_1)
if not options.dir_fastq_2:
	sys.exit("Missing input fastq 2 file, -2 xxx.fastq")
if not os.path.isfile(options.dir_fastq_1):
	sys.exit("Missing input fastq 2 file, %r" % options.dir_fastq_2)
if not options.dir_out:
	sys.exit("Missing output directory, -o xxx.fasta")

dir_fastq_1 = os.path.abspath(options.dir_fastq_1)
dir_fastq_2 = os.path.abspath(options.dir_fastq_2)
dir_out_fa = options.dir_out

print("begin convert")

if dir_fastq_1.endswith(".gz"):
	with gzip.open(dir_fastq_1, 'rb') as f:
		fastq_1 = f.read().split('\n')
else:
	fastq_1 = open(dir_fastq_1, 'r').read().split('\n')

if dir_fastq_2.endswith(".gz"):
	with gzip.open(dir_fastq_2, 'rb') as f:
		fastq_2 = f.read().split('\n')
else:
	fastq_2 = open(dir_fastq_2, 'r').read().split('\n')

fastq_1 = fastq_1[:-1]
fastq_2 = fastq_2[:-1]

print("length of fastq1: " + str(len(fastq_1)))
print("length of fastq2: " + str(len(fastq_2)))

w_file = open(dir_out_fa, 'w')
for i in range(0,len(fastq_1),4):
	w_file.write(">" + fastq_1[i][1:-8] + ".1\n" + fastq_1[i+1] + "\n")
	w_file.write(">" + fastq_2[i][1:-8] + ".2\n" + fastq_2[i+1] + "\n")
