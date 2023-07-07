#Get regions where there are reads < k bp
import pysam
import sys 
import os
import numpy as np
import pybedtools

#TODO: optimize SV filters

#Get regions on ref where its not covered by at least one of the assembly
def get_non_cover_regions(assem_bam_file):
    samfile = pysam.AlignmentFile(assem_bam_file, "rb")
    chr_list = samfile.header.references
    for chr_name in chr_list:
        #test
        #print(chr_name)
        #stop when parsed no_of_reads
        #no_of_reads = 400000
        cur_end = 0
        #loop through iter, index will not be reset
        #if if_hg38:
        #    ref_name = "chr" + chr_name
        #else:
        #    ref_name = chr_name
        ref_name = chr_name
        iter = samfile.fetch(ref_name)

        for rec in iter:
            if rec.reference_start > cur_end:
                print(str(ref_name) + "\t" + str(cur_end) + "\t" + str(rec.reference_start))
                cur_end = rec.reference_end
            else:
                if rec.reference_end > cur_end:
                    cur_end = rec.reference_end
    
#Get SV positions
# TODO update this too (LD)
def get_sv_positions(output_dir, vcf_file):
    f = pysam.VariantFile(vcf_file, 'r')
    #create bedfile contains SVs' positions
    g = open(output_dir + "SV_positions.bed", "w")
    for counter, rec in enumerate(f.fetch()):
        ref_name = rec.chrom
        sv_type = rec.info['SVTYPE']
        sv_len = rec.rlen
        #TODOL double check the start for different types
        sv_pos = rec.pos
        sv_end = rec.stop
        if sv_type not in ['DEL', 'INS', 'INV', 'DUP']:
            continue
        g.write(str(ref_name) + "\t")
        g.write(str(sv_pos) + "\t")
        g.write(str(sv_end))
        g.write("\n")
    g.close()
    f.close()

#Output filtered calls' info in non-covered regions
def output_non_cov_call_info(output_dir, SV_positions_file, assem1_non_cov_regions_file, assem2_non_cov_regions_file):
    SV_positions = pybedtools.BedTool(SV_positions_file)
    assem1_non_cov_regions = pybedtools.BedTool(assem1_non_cov_regions_file)
    assem2_non_cov_regions = pybedtools.BedTool(assem2_non_cov_regions_file)
    
    exclude_assem1_non_cover = SV_positions.intersect(assem1_non_cov_regions, u = True)
    exclude_assem2_non_cover = SV_positions.intersect(assem2_non_cov_regions, u = True)
    
    exclude_assem1_non_cover.saveas(output_dir + 'exclude_assem1_non_cover.bed')
    exclude_assem2_non_cover.saveas(output_dir + 'exclude_assem2_non_cover.bed')
    

#Get regions where read depth > 2 * avg_read_depth
#For now, we filter calls by read depth
#In other words, here output calls having high read depth
#TODO: Too slow here
def get_high_depth_calls_info(output_dir, read_bam_file, vcf_file, avg_read_depth):
    sv_len_limit = 100000
    avg_read_depth = float(avg_read_depth)
    
    samfile = pysam.AlignmentFile(read_bam_file, "rb")
    f = pysam.VariantFile(vcf_file,'r')
    #TODO: change open condition
    g = open(output_dir + "exclude_high_depth.bed", "w")
    for counter, rec in enumerate(f.fetch()):
        #test
        #if counter > 250:
        #    break
        
        #get ref start and ref end
        name = rec.chrom
        sv_pos = rec.pos
        sv_end = rec.stop
        sv_type = rec.info['SVTYPE']
        if sv_type not in ['DEL', 'INS', 'INV', 'DUP']:
            continue
        
        sv_len = abs(rec.info['SVLEN'][0]) 
        if sv_len > sv_len_limit:
            continue
        
        #print(sv_end - sv_pos)
        res = samfile.count_coverage(name, sv_pos, sv_end+1, quality_threshold = 0)
        #print(res)
        if round(np.sum(res)/(sv_end+1-sv_pos), 2) > 2*avg_read_depth:
            #test
            #print(counter, round(np.sum(res)/(sv_end+1-sv_pos), 2))
            g.write(str(name) + "\t")
            g.write(str(sv_pos) + "\t")
            g.write(str(sv_end))
            g.write("\n")
    g.close()
    f.close()

def main():
    #get command line input
    #n = len(sys.argv)
    #assembly bam file
    bam_file = sys.argv[1]
#     #if_hg38 = True
#     avg_read_depth = sys.argv[5]
#     #reads bam file
#     read_bam_file = sys.argv[6]
#     #callset file
#     vcf_file = sys.argv[7]

    #Output regions on ref where its not covered by at least one of the assembly
    get_non_cover_regions(bam_file)
    
    #Output sv positions
#     get_sv_positions(output_dir, vcf_file)
    
    #Output filtered calls in non-covered regions
#     SV_positions_file = output_dir + "SV_positions.bed"
#     assem1_non_cov_regions_file = output_dir + "assem1_non_cov_regions.bed"
#     assem2_non_cov_regions_file = output_dir + "assem2_non_cov_regions.bed"
#     output_non_cov_call_info(output_dir, SV_positions_file, assem1_non_cov_regions_file, assem2_non_cov_regions_file)
    
    #Get regions where read depth > 2 * avg_read_depth
#     get_high_depth_calls_info(output_dir, read_bam_file, vcf_file, avg_read_depth)
    
if __name__ == "__main__":
    main()