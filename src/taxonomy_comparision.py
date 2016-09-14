import sys, os, stat
from optparse import OptionParser

# # Parse Command Line
# usage = """
# Description: taxonomy comparsion tool 

# """

# parser = OptionParser(usage=usage)
# parser.add_option("-i", "--input", dest="dir_input",
# 	default=None, metavar="FILE",
# 	help="input directory (required)")
# parser.add_option("-o", "--output", dest="dir_out",
# 	default=None, metavar="FILE",
# 	help="Output filename (required)")

# options,args = parser.parse_args()

# if not options.dir_input:
# 	sys.exit("Missing input fasta file, -i dir")
# if not os.path.isfile(options.dir_input):
# 	sys.exit("Missing input fasta file, %r" % options.dir_input)
# if not options.dir_out:
# 	sys.exit("Missing output directory, -o xxxxx")
# if not options.node:
# 	sys.exit("Missing number of node")
# list_input_file = ["visualize_taxonomy/GOEL01S1L001.level_7.converted.txt", "visualize_taxonomy/GOEL02S2L001.level_7.converted.txt" ]
list_input_file = sys.argv[1:]
short_name = []
sum_species = []
for file_name in list_input_file:
	file_name = file_name.split("/")[-1].split(".")[0]
	short_name.append(file_name)


list_taxa = {}
for i in range(len(list_input_file)):
	input_file = os.path.abspath(list_input_file[i])
	num_species = 0
	file = open(input_file, 'r')
	for line in file:
		num = int(line[:line.find('\t')])
		taxa = line[line.find('\t')+1:-1]
		if taxa not in list_taxa.keys():
			list_taxa[taxa]=[]
			for j in range(len(list_input_file)):
				list_taxa[taxa].append(0)
			list_taxa[taxa][i] = num
			num_species+=num
		else:
			list_taxa[taxa][i] = num
			num_species+=num
	sum_species.append(num_species)

for name in short_name:
	sys.stdout.write(name + "\t")
sys.stdout.write("taxonomic\n")

for taxa, num in list_taxa.items():
	for num_taxa in num:
		sys.stdout.write(str(num_taxa) + "\t")
	sys.stdout.write(taxa.replace("\t","__") + "\n")

for num_species in sum_species:
	sys.stdout.write(str(num_species) + "\t")
sys.stdout.write("\n")
