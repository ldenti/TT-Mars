# %%
#generate bed file that contains positions on ref to be lifted to contigs
import csv
import math
#pysam: https://pysam.readthedocs.io/en/latest/usage.html
import pysam
#print(pysam.__version__)
import sys 
  
#get command line input
#n = len(sys.argv)
assembly_bam_file_hap1 = sys.argv[1]
#print("\nName of Python script:", sys.argv[0])

interval = 20

#check if k-th bit of a given number is set or not 
def isKthBitSet(n, k): 
    if n & (1 << (k - 1)): 
        return True 
    else: 
        return False

#for assem1
#generate bed file liftover

samfile = pysam.AlignmentFile(assembly_bam_file_hap1, "rb")
for record in samfile:
    #TODO: solve the unmapped segment problem???
    if not isKthBitSet(record.flag, 3):
        #samLiftover's bed file needs ref name for both directions
        ref_name = record.reference_name
        query_name = record.query_name
        #TODO: check not matching with paf
        query_start = record.query_alignment_start
        query_end = record.query_alignment_end
        ref_start = record.reference_start
        ref_end = record.reference_end
        query_start = math.ceil(query_start/interval) * interval
        query_end = math.floor(query_end/interval) * interval
        ref_start = math.ceil(ref_start/interval) * interval
        ref_end = math.floor(ref_end/interval) * interval
        
        #if start with hard clip
        if record.cigartuples[0][0] == 5:
            query_start += record.cigartuples[0][1]
            query_end += record.cigartuples[0][1]

        if query_end < query_start:
            tem = query_start
            query_start = query_end
            query_end = tem

        for i in range(query_start, query_end + interval, interval):
            print(query_name + "\t" + str(i) + "\t" + str(i+1) + "\t" + ref_name)
