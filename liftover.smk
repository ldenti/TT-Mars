from os.path import join as pjoin

REF = config["fa"]
HAP1 = config["hap1"]
HAP2 = config["hap2"]
ODIR = config["odir"]

rule run:
    input:
        expand(pjoin(ODIR, "hap{H}.non_covered_regions.bed"),
            H = ["1", "2"]),
        expand(pjoin(ODIR, "hap{H}.lo_pos_result.{_}.compressed.bed"),
            H = ["1", "2"],
            _ = ["0", "1"]),

rule index_ref:
    input:
        fa = REF
    output:
        mms=REF + ".mms",
        gli=REF + ".gli",
    # conda: "envs/lra.yaml"
    threads: 1
    shell:
        """
        lra index -CONTIG {input.fa}
        """

rule align:
    input:
        fa=REF,
        mms=REF + ".mms",
        gli=REF + ".gli",
        hap=lambda wildcards: HAP1 if wildcards.H == "1" else HAP2
    output:
        bam=pjoin(ODIR, "hap{H}.bam")
    # conda: "envs/lra.yaml"
    threads: 16
    shell:
        """
        lra align -CONTIG {input.fa} {input.hap} -t {threads} -p s | samtools sort > {output.bam}
        samtools index {output.bam}
        """

rule trim:
    input:
        bam = rules.align.output.bam
    output:
        bam=pjoin(ODIR, "hap{H}.nooverlap.bam")
    threads: 1
    shell:
        """
        python3 trim_overlapping_contigs.py {input.bam} | samtools sort > {output.bam}
        samtools index {output.bam}
        """

rule bam2sam:
    input:
        bam=rules.trim.output.bam,
    output:
        sam=pjoin(ODIR, "hap{H}.nooverlap.sam")
    shell:
        """
        samtools view -h {input.bam} > {output.sam}
        """

rule bam2bed:
    input:
        bam=rules.trim.output.bam,
    output:
        bed=pjoin(ODIR, "hap{H}.lo_pos.1.bed"),
    shell:
        """
        python3 lo_assem_to_ref.py {input.bam} > {output.bed}
        """

rule samliftover:
    input:
        sam=rules.bam2sam.output.sam,
        bed=rules.bam2bed.output.bed,
    output:
        bed=pjoin(ODIR, "hap{H}.lo_pos_result.1.bed")
    shell:
        """
        /home/denti/software/mcutils/src/samLiftover {input.sam} {input.bed} {output.bed} --dir 1
        """

rule bam2bed_0:
    input:
        bam=rules.trim.output.bam,
    output:
        bed=pjoin(ODIR, "hap{H}.lo_pos.0.bed"),
    shell:
        """
        python3 lo_assem_to_ref_0.py {input.bam} > {output.bed}
        """

rule samliftover_0:
    input:
        sam=rules.bam2sam.output.sam,
        bed=rules.bam2bed_0.output.bed,
    output:
        bed=pjoin(ODIR, "hap{H}.lo_pos_result.0.bed")
    shell:
        """
        /home/denti/software/mcutils/src/samLiftover {input.sam} {input.bed} {output.bed} --dir 0
        """

rule compress_liftover:
    input:
        bed=pjoin(ODIR, "hap{H}.lo_pos_result.{_}.bed")
    output:
        bed=pjoin(ODIR, "hap{H}.lo_pos_result.{_}.compressed.bed")
    shell:
        """
        python3 compress_liftover.py {input.bed} > {output.bed}
        """

rule get_noncovered_regions:
    input:
        bam=rules.trim.output.bam,
    output:
        bed=pjoin(ODIR, "hap{H}.non_covered_regions.bed")
    shell:
        """
        python3 get_conf_int.py {input.bam} > {output.bed}
        """