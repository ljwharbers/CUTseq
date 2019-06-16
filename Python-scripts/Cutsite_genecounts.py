"""

Get the amount of cuts per gene

"""

from Bio import SeqIO
from Bio.Restriction import *

#define variables
#enzyme = {"BfaI", "CviAII", "FatI", "NlaIII", "MseI", "AflII", "AseI", "BsrGI", "HindIII", "NcoI", "NdeI", "NsiI",
#          "PciI", "PstI", "SacI", "SpeI", "SphI", "EcoRI"}

commercial = CommOnly

enzyme = []
for i in commercial:
    if i.is_blunt() == False:
        if len(i) == 4:             #change len(i) to len(site.i) ??
            enzyme.append(i)
reference = "E:\SciLife Lab\CUTseq\Data\References/cosmic_exons_padded_wochr.test.fa"

#run either of the outputs depending on whats needed

output = open("E:\SciLife Lab\CUTseq\Cut_site_analysis/rb_batch_Subset_cutsites_newfasta.txt", "w")

rb = RestrictionBatch(enzyme)

for record in SeqIO.parse(reference, "fasta"):
    sites = rb.search(record.seq)               #search for cut site overlaps
    for x in sites:                             #get overlap and select the regions with len() > 0
        if(len(sites[x])) > 0:
           # write record ID, region and amount of cut sites
           output.write(record.description + "\t" + str(x) + "\t" + str(len(sites[x])) + "\n")