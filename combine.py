import csv
import sys
import os

import argparse

import fn
import ttmars

# parser = argparse.ArgumentParser()
# parser.add_argument("output_dir",
#                     help="output directory")
# parser.add_argument("no_X_chr",
#                     choices=[1, 2],
#                     help="male sample 1, female sample 2",
#                     type=int)
# parser.add_argument("-v",
#                     "--vcf_out",
#                     help="output results as vcf files, must be used together with -f/--vcf_file",
#                     action="store_true")
# parser.add_argument("-f",
#                     "--vcf_file",
#                     help="input vcf file using as template, must be used together with -v/--vcf_out or -n/--false_neg")
# parser.add_argument("-g",
#                     "--gt_vali",
#                     help="conduct genotype validation",
#                     action="store_true")
# parser.add_argument("-n",
#                     "--false_neg",
#                     help="output false negative, must be used together with -t/--truth_file and -f/--vcf_file",
#                     action="store_true")
# parser.add_argument("-t",
#                     "--truth_file",
#                     help="input truth vcf file, must be used together with -n/--false_neg")
# args = parser.parse_args()

# if bool(args.vcf_out) ^ bool(args.vcf_file):
#     parser.error('-f/--vcf_file must be used with -v/--vcf_out')
    
# if bool(args.false_neg) ^ bool(args.truth_file):
#     parser.error('-t/--truth_file must be used with -n/--false_neg')
    
# if (bool(args.false_neg)) and (not bool(args.vcf_file)):
#     parser.error('-n/--false_neg must be used with -f/--vcf_file')

# output_dir = args.output_dir + "/"

# if int(args.no_X_chr) == 1:
#     if_male = True
# else:
#     if_male = False
    
# if args.vcf_out:
#     if_vcf = True
#     in_vcf_file = args.vcf_file
# else:
#     if_vcf = False
    
# if args.false_neg:
#     output_fn = True
#     in_truth_file = args.truth_file
# else:
#     output_fn = False

output_dir = ttmars.output_dir
if_male = ttmars.if_male
if_vcf = ttmars.if_vcf
in_vcf_file = ttmars.vcf_file
output_fn = ttmars.output_fn
if output_fn:
    in_truth_file = ttmars.in_truth_file
if_gt = ttmars.if_gt

def combine_output():
    other_sv_res_file = output_dir+"ttmars_res.txt"
    regdup_res_file = output_dir+"ttmars_regdup_res.txt"
    agg_res_file = output_dir+"ttmars_agg_res.txt"

    with open(other_sv_res_file) as f:
        reader = csv.reader(f, delimiter="\t")
        other_sv_res = list(reader)
    f.close()

    if_dup = True
    try:
        with open(regdup_res_file) as f:
            reader = csv.reader(f, delimiter="\t")
            regdup_res = list(reader)
        f.close()
    except:
        if_dup = False

    if_agg = True
    try:
        with open(agg_res_file) as f:
            reader = csv.reader(f, delimiter="\t")
            agg_res = list(reader)
        f.close()
    except:
        if_agg = False


    #if output fn
    if output_fn:
        #input truth set
        sv_list = fn.idx_sv(in_truth_file)
        #input candidate set
        cand_sv_list = fn.idx_sv(in_vcf_file)

        print(len(sv_list), len(cand_sv_list))

        sv_list_sorted = fn.sort_sv_list(sv_list)
        cand_sv_list_sorted = fn.sort_sv_list(cand_sv_list)

        tp_base_ctr = fn.count_tp_base_dist_only(sv_list_sorted, cand_sv_list_sorted)
        recall = tp_base_ctr / len(sv_list_sorted)

        print(tp_base_ctr, len(sv_list_sorted))
        print("Recall of candidate callset: " + str(recall))

    sv_dict = {}
    #results of SVs other than interspersed DUP
    for rec in other_sv_res:
        sv_idx = rec[0]
        ref_name = rec[1]
        sv_pos = int(rec[2])
        sv_end = int(rec[3])
        sv_type = rec[4]

        if not if_gt:
            sv_dict[sv_idx] = [rec[5], rec[6], rec[7], ref_name, int(sv_pos), int(sv_end), sv_type]
        else:
            sv_dict[sv_idx] = [rec[5], rec[6], rec[7], ref_name, int(sv_pos), int(sv_end), sv_type, rec[8]]

    #     if not if_gt:
    #         sv_dict[(ref_name, int(sv_pos), int(sv_end), sv_type)] = [rec[4], rec[5], rec[6]]
    #     else:
    #         sv_dict[(ref_name, int(sv_pos), int(sv_end), sv_type)] = [rec[4], rec[5], rec[6], rec[7]]

    #interspersed DUP
    if if_dup:
        for rec in regdup_res:
            sv_idx = rec[0]
            ref_name = rec[1]
            sv_pos = int(rec[2])
            sv_end = int(rec[3])
            sv_type = rec[4]

            if sv_idx in sv_dict:
                if not if_gt:
                    if rec[7] == 'True' and sv_dict[sv_idx][2] == 'False':
                        sv_dict[sv_idx] = [rec[5], rec[6], rec[7], ref_name, int(sv_pos), int(sv_end), sv_type]
                else:
                    if rec[7] == 'True' and sv_dict[sv_idx][2] == 'False':
                        sv_dict[sv_idx] = [rec[5], rec[6], rec[7], ref_name, int(sv_pos), int(sv_end), sv_type, rec[8]]
                    elif rec[7] == 'True' and sv_dict[sv_idx][2] == 'True':
                        if rec[8] == 'True' and sv_dict[sv_idx][3] == 'False':
                            sv_dict[sv_idx] = [rec[5], rec[6], rec[7], ref_name, int(sv_pos), int(sv_end), sv_type, rec[8]]
            else:
                if not if_gt:
                    sv_dict[sv_idx] = [rec[5], rec[6], rec[7], ref_name, int(sv_pos), int(sv_end), sv_type]
                else:
                    sv_dict[sv_idx] = [rec[5], rec[6], rec[7], ref_name, int(sv_pos), int(sv_end), sv_type, rec[8]]

    #agg sv
    if if_agg:
        for rec in agg_res:
            sv_idx = rec[0]
            ref_name = rec[1]
            sv_pos = int(rec[2])
            sv_end = int(rec[3])
            sv_type = rec[4]

            if sv_idx in sv_dict:
                if not if_gt:
                    if rec[7] == 'True' and sv_dict[sv_idx][2] == 'False':
                        sv_dict[sv_idx] = [rec[5], rec[6], rec[7], ref_name, int(sv_pos), int(sv_end), sv_type]
                else:
                    if rec[7] == 'True' and sv_dict[sv_idx][2] == 'False':
                        sv_dict[sv_idx] = [rec[5], rec[6], rec[7], ref_name, int(sv_pos), int(sv_end), sv_type, rec[8]]
                    elif rec[7] == 'True' and sv_dict[sv_idx][2] == 'True':
                        if rec[8] == 'True' and sv_dict[sv_idx][3] == 'False':
                            sv_dict[sv_idx] = [rec[5], rec[6], rec[7], ref_name, int(sv_pos), int(sv_end), sv_type, rec[8]]
            else:
                if not if_gt:
                    if rec[7] == 'True':
                        sv_dict[sv_idx] = [rec[5], rec[6], rec[7], ref_name, int(sv_pos), int(sv_end), sv_type]
                else:
                    if rec[7] == 'True':
                        sv_dict[sv_idx] = [rec[5], rec[6], rec[7], ref_name, int(sv_pos), int(sv_end), sv_type, rec[8]]

    if if_male:
        chrx_res_file = output_dir+"ttmars_chrx_res.txt"
        chrx_res = []
        if os.path.exists(chrx_res_file):
            with open(chrx_res_file) as f:
                reader = csv.reader(f, delimiter="\t")
                chrx_res = list(reader)
        f.close()
        for rec in chrx_res:
            if len(rec) == 0:
                continue

            sv_idx = rec[0]
            ref_name = rec[1]
            sv_pos = int(rec[2])
            sv_end = int(rec[3])
            sv_type = rec[4]

            if sv_idx in sv_dict:
                if not if_gt:
                    if rec[7] == 'True' and sv_dict[sv_idx][2] == 'False':
                        sv_dict[sv_idx] = [rec[5], rec[6], rec[7], ref_name, int(sv_pos), int(sv_end), sv_type]
                else:
                    if rec[7] == 'True' and sv_dict[sv_idx][2] == 'False':
                        sv_dict[sv_idx] = [rec[5], rec[6], rec[7], ref_name, int(sv_pos), int(sv_end), sv_type, rec[8]]
                    elif rec[7] == 'True' and sv_dict[sv_idx][2] == 'True':
                        if rec[8] == 'True' and sv_dict[sv_idx][3] == 'False':
                            sv_dict[sv_idx] = [rec[5], rec[6], rec[7], ref_name, int(sv_pos), int(sv_end), sv_type, rec[8]]
            else:
                if not if_gt:
                    sv_dict[sv_idx] = [rec[5], rec[6], rec[7], ref_name, int(sv_pos), int(sv_end), sv_type]
                else:
                    sv_dict[sv_idx] = [rec[5], rec[6], rec[7], ref_name, int(sv_pos), int(sv_end), sv_type, rec[8]]  

    g = open(output_dir + "/ttmars_combined_res.txt", "w")
    for key in sv_dict:
        res = sv_dict[key]

        g.write(str(key) + "\t")
        g.write(str(res[0]) + "\t")
        g.write(str(res[1]) + "\t")
        g.write(str(res[2]) + "\t")
        g.write(str(res[3]) + "\t")
        g.write(str(res[4]) + "\t")
        g.write(str(res[5]) + "\t")
        g.write(str(res[6]))

        if if_gt:
            g.write("\t" + str(res[7]))

        g.write("\n")
    g.close()

    #if output vcf
    if if_vcf:
        from pysam import VariantFile
        vcf_in = VariantFile(in_vcf_file)
        vcfh = vcf_in.header
        #vcfh.add_meta('INFO', items=[('ID',"TTMars"), ('Number',1), ('Type','String'),('Description','TT-Mars NA12878 results: TP, FA, NA or .')])
        vcfh.add_meta('INFO', items=[('ID',"GT_vali"), ('Number',1), ('Type','String'),('Description','TT-Mars GT validation (require flag -g): True, False or NA')])
        vcf_out_tp = VariantFile(output_dir+"ttmars_tp.vcf", 'w', header=vcfh)
        vcf_out_fp = VariantFile(output_dir+"ttmars_fp.vcf", 'w', header=vcfh)
        vcf_out_na = VariantFile(output_dir+"ttmars_na.vcf", 'w', header=vcfh)

        for count, rec in enumerate(vcf_in.fetch()):
            key = str(count)
            try:
                validation_res = sv_dict[key][2]
                if validation_res == 'True':
                    if if_gt:
                        rec.info['GT_vali'] = sv_dict[key][7]
                    vcf_out_tp.write(rec)
                elif validation_res == 'False':
                    vcf_out_fp.write(rec)
            except:
                vcf_out_na.write(rec)

    #     for rec in vcf_in.fetch():
    #         ref_name = rec.chrom
    #         sv_type = rec.info['SVTYPE']
    #         sv_pos = rec.pos
    #         sv_end = rec.stop

    #         try:
    #             validation_res = sv_dict[(ref_name, int(sv_pos), int(sv_end), sv_type)][2]
    #             if validation_res == 'True':
    #                 if if_gt:
    #                     rec.info['GT_vali'] = sv_dict[(ref_name, int(sv_pos), int(sv_end), sv_type)][3]
    #                 vcf_out_tp.write(rec)
    #             elif validation_res == 'False':
    #                 vcf_out_fp.write(rec)
    #         except:
    #             vcf_out_na.write(rec)

        vcf_out_tp.close()
        vcf_out_fp.close()
        vcf_out_na.close()
        vcf_in.close()    
    
#remove files
def remove_files():
    import os
    for name in ['assem1_non_cov_regions.bed', 'assem2_non_cov_regions.bed',
                 'exclude_assem1_non_cover.bed', 'exclude_assem2_non_cover.bed',
                 'SV_positions.bed', 'ttmars_chrx_res.txt', 'ttmars_regdup_res.txt',
                 'ttmars_res.txt', 'align_info_assem1_chrall.txt', 'align_info_assem2_chrall.txt',
                 'all_reg_dup.fasta', 'all_reg_dup.fasta.fai']:
        if os.path.exists(output_dir + name):
            os.remove(output_dir + name)








