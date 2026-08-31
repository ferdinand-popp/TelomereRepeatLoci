# Author: Lina Sieverling

#!/usr/bin/python

# Usage: source activate telomereEnv
#        python find_discordant_reads.py \
#                -i <*_filtered_intratelomeric.bam resulting from TelomereHunter> \
#                -o <*_tumor_telomere_insertions.tsv output file>
# Description: extracts read names and positions of mates for reads that fullfill the following criteria:
#               - 1 mate is intratelomeric, the other mate is not
#               - mate that is not intratelomeric needs to be mapped (attention: no information about mapping quality)



import os
import sys, getopt
import re
import pysam

# ----------------------------------------------------------------
# read command line args 
# ----------------------------------------------------------------
myopts, args = getopt.getopt(sys.argv[1:],"i:o:")
 
for opt, arg in myopts:
  if opt == '-i':
    intratel_bam = arg
  elif opt == '-o':
    outfile_path = arg
  else:
    print("Usage: %s -i input_path (TelomereHunter intratelomeric bam file) -o outfile_path" % sys.argv[0])

#####################################################################################################################################


print("[find_discordant_reads] input bam: %s" % intratel_bam)
print("[find_discordant_reads] output tsv: %s" % outfile_path)

# ----------------------------------------------------------------
# make a dictionary with the number of reads with the same name
# ----------------------------------------------------------------
read_name_dict = {}

bamfile = pysam.Samfile( intratel_bam, "rb" )

total_reads = 0
for read in bamfile.fetch(until_eof=True):

  total_reads += 1
  read_name = read.qname

  try:
    read_name_dict[read_name] += 1
  except:
    read_name_dict[read_name] = 1

bamfile.close()

print("[find_discordant_reads] scanned %d reads in intratelomeric bam (%d distinct read names)" % (total_reads, len(read_name_dict)))



# ---------------------------------------------------------------------------
# go through bam file again and extract mate mapping positions of reads if
# the mate is mapped and not intratelomeric
# ---------------------------------------------------------------------------
bamfile2 = pysam.Samfile( intratel_bam, "rb" )

output = "read_name\tmate_chr\tmate_position\n"

n_skipped_mate_unmapped = 0
n_skipped_mate_ref_unknown = 0
n_skipped_mate_intratelomeric = 0
n_discordant = 0

for read in bamfile2.fetch(until_eof=True):

  #skip reads where the mate is unmapped
  if read.mate_is_unmapped:
    n_skipped_mate_unmapped += 1
    continue

  #skip reads where the reference ID of the mate is not known ('*', this can happen when the mapq of the mate is 255='not known')
  if read.next_reference_id==-1:
    n_skipped_mate_ref_unknown += 1
    continue

  #skip reads where the mate is also intratelomeric
  read_name = read.qname
  if read_name_dict[read_name] == 2:
    n_skipped_mate_intratelomeric += 1
    continue

  #get chromosome of mate
  mate_chr = read.next_reference_name

  # get 0-based mapping position of mate, adding 1 to get it 1-based like in SAM file
  mate_position = read.next_reference_start + 1

  output += read_name + "\t" + str(mate_chr) + "\t" + str(mate_position) + "\n"
  n_discordant += 1

bamfile2.close()

print("[find_discordant_reads] skipped %d reads with unmapped mate" % n_skipped_mate_unmapped)
print("[find_discordant_reads] skipped %d reads with unknown mate reference" % n_skipped_mate_ref_unknown)
print("[find_discordant_reads] skipped %d reads with intratelomeric mate" % n_skipped_mate_intratelomeric)
print("[find_discordant_reads] found %d discordant reads" % n_discordant)


# ----------------------------------------------------------------
# write read name and mapping position of mate to output table
# ----------------------------------------------------------------
outfile = open( outfile_path, "w")
outfile.write(output)
outfile.close()

print("[find_discordant_reads] wrote %d discordant reads to %s" % (n_discordant, outfile_path))





