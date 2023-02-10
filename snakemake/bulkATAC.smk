configfile: "config.yaml"

# /bin/nice -n5 snakemake -s /home/dmakosa/working_data_04/Bulk_ATAC_QCmetrics_tool/snakemake/bulkATAC.smk --use-conda --default-resources "tmpdir='/scratchfs/dmakosa/tmp'" --cores 64
# snakemake --forceall --rulegraph -s test2.smk | dot -Tpdf > dag.pdf

# ---- DICTIONARIES ---- #
READS = ["R1", "R2"]
LOGs = ["logPreQCs", "logAdap", "logPostQCs"]
DIRS = ["output/1.QCs", "output/3.QCsRmAdaptors/"]
histogram_script = "/home/dmakosa/working_data_04/Bulk_ATAC_QCmetrics_tool/dependencies/Frag_hist.R"

index = config["index"] 
blacklistedMito = config["blacklistedMito"] 
blacklisted = config["blacklisted"] 
promoters = config["promoters"] 
size = config["size"] 
generalPeakLanscape = config["generalPeakLanscape"]
OliviersPeakLanscape = config["OliviersPeakLanscape"] 


# ---- TARGET RULE ---- # 
rule all:
    input:
        "output/5.SubsamplingOfFiltered/4.summary/summary.tsv",
        expand("output/5.SubsamplingOfFiltered/3.readsInPeaks/NumberOfReads_{sample}.txt", sample=config["samples"]),
        expand("output/5.SubsamplingOfFiltered/3.readsInPeaks/ReadsInPeaks_{sample}.txt", sample=config["samples"]),
        expand("output/5.SubsamplingOfFiltered/2.callpeaks/NFRPeaks_{sample}.txt", sample=config["samples"]),
        expand("output/5.SubsamplingOfFiltered/1.subsampling/subsampled_{sample}.txt", sample=config["samples"]),
        expand("output/5.SubsamplingOfFiltered/0.proportions/{sample}.proportions", sample=config["samples"]),
        expand("output/4.Alignment/flagstatsCombined/{sample}.flagStats_combined", sample=config["samples"]),
        expand("output/4.Alignment/mapped_pairs_dedup_filt_noMT/{sample}.noOfReads", sample=config["samples"]),
        expand("output/4.Alignment/readsInProm/{sample}_pairs_dedup_filt_noMT.FRIPromoters", sample=config["samples"]),
        expand("output/4.Alignment/OliviersPeakLanscape/{sample}_pairs_dedup_filt_noMT.FRIPOlivier", sample=config["samples"]),
        expand("output/4.Alignment/OliviersPeakLanscape/{sample}_pairs_dedup_filt_noMT.ReadsInOliviersPeakLanscape", sample=config["samples"]),
        expand("output/4.Alignment/OliviersPeakLanscape/{sample}_pairs_dedup_filt_noMT.ReadsInOliviersPeakLanscape_overProm", sample=config["samples"]),
        expand("output/4.Alignment/GeneralPeakLandscape/{sample}_pairs_dedup_filt_noMT.FRIP", sample=config["samples"]),
        expand("output/4.Alignment/InsertSizes/{sample}_pairs_dedup_filt_noMT.InsertSizesBAMPE_hist.png", sample=config["samples"]),
        expand("output/4.Alignment/InsertSizes/{sample}_pairs_dedup_filt_noMT.InsertSizesBAMPE", sample=config["samples"]),
        expand("output/4.Alignment/Peaks/{sample}_pairs_dedup_filt_noMT.ReadsInNFRpeaks_overProm", sample=config["samples"]),
        expand("output/4.Alignment/GeneralPeakLandscape/{sample}_pairs_dedup_filt_noMT.ReadsInGeneralPeakLandscape", sample=config["samples"]),
        expand("output/4.Alignment/GeneralPeakLandscape/{sample}_pairs_dedup_filt_noMT.ReadsInGeneralPeakLandscape_overProm", sample=config["samples"]),
        expand("output/4.Alignment/Peaks/{sample}_pairs_dedup_filt_noMT.ReadsInBAMPEpeaks", sample=config["samples"]),
        expand("output/4.Alignment/Peaks/{sample}_pairs_dedup_filt_noMT.ReadsInNFRpeaks_overProm", sample=config["samples"]),
        expand("output/4.Alignment/Peaks/{sample}_pairs_dedup_filt_noMT.ReadsInNFRpeaks", sample=config["samples"]),
        expand("output/4.Alignment/Peaks/{sample}_pairs_dedup_filt_noMT.BAMPE_peaks.narrowPeak", sample=config["samples"]),
        expand("output/4.Alignment/Peaks/{sample}_pairs_dedup_filt_noMT.NFR_peaks.narrowPeak", sample=config["samples"]),
        expand("output/4.Alignment/readsInProm/{sample}_pairs_dedup_filt_noMT.ReadsInProm", sample=config["samples"]),
        expand("output/4.Alignment/coverage/{sample}_pairs_dedup_filt_noMT.coverage.bigwig", sample=config["samples"]),
        expand("output/4.Alignment/mapped_pairs_dedup_filt_mitocounts/{sample}_pairs_dedup_filt.mitoReads", sample=config["samples"]),
        expand("output/4.Alignment/mapped_pairs_dedup_filt_noMT/{sample}_pairs_dedup_filt_noMT.flagStats", sample=config["samples"]),
        expand("output/4.Alignment/mapped_pairs_dedup_filt/{sample}_pairs_dedup_filt.flagStats", sample=config["samples"]),
        expand("output/4.Alignment/mapped_pairs_dedup/{sample}_pairs_dedup.flagStats", sample=config["samples"]),
        expand("output/4.Alignment/mapped_pairs/{sample}_pairs.flagStats", sample=config["samples"]),
        expand("output/4.Alignment/mapped/{sample}.flagStats", sample=config["samples"]),
        expand("output/4.Alignment/mapped/{sample}.bam", sample=config["samples"]),
        expand("output/2.ProcessedReads/{sample}_R1_clean.fastq.gz", sample=config["samples"]),
        expand("output/2.ProcessedReads/{sample}_R2_clean.fastq.gz", sample=config["samples"]),
        expand("output/0.LOGs/logPreQCs_{sample}.log", sample=config["samples"]),
        expand("output/0.LOGs/logPostQCs_{sample}.log", sample=config["samples"]),
        expand("{dir}", dir=DIRS),
        expand("data/{sample}_{read}_001.fastq.gz", sample=config["samples"], read=READS)
    message: "Target rule"



# ---- RULES ---- # 

# ---- PRE-PROCESSING ---- # 
rule makedir:
    output: 
        dir1 = directory("output/1.QCs"),
        dir3 = directory("output/3.QCsRmAdaptors/")
    message: "creating directories..."
    threads: 1
    shell: "mkdir -p {output.dir1} {output.dir3}"
    
rule fastqc_raw_reads:
    input:
        fq1 = "data/{sample}_R1_001.fastq.gz",
        fq2 = "data/{sample}_R2_001.fastq.gz",
        directory_check = "output/1.QCs",
        directory_check2 = "output/3.QCsRmAdaptors"
    output: "output/0.LOGs/logPreQCs_{sample}.log"
    message: "QC of raw reads on sample {wildcards.sample} files {input.fq1} {input.fq2}"
    threads: 1
    conda:
        "10xmethylomes"
    shell: "fastqc --quiet -o output/1.QCs/ {input.fq1} {input.fq2} > {output}"

rule clean_adapters:
    input:
        fq1 = "data/{sample}_R1_001.fastq.gz",
        fq2 = "data/{sample}_R2_001.fastq.gz",
        directory_check = "output/1.QCs",
        directory_check2 = "output/3.QCsRmAdaptors"
    output:
        fq1 = "output/2.ProcessedReads/{sample}_R1_clean.fastq.gz",
        fq2 = "output/2.ProcessedReads/{sample}_R2_clean.fastq.gz"
    log: "output/0.LOGs/logAdap_{sample}.log"
    message: "Trimming adapters on sample {wildcards.sample} files {input.fq1} {input.fq2}"
    threads: 1
    conda:
        "10xmethylomes"
    shell:
        "cutadapt -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -A CTGTCTCTTATACACATCTGACGCTGCCGACGA "
        "-o {output.fq1} -p {output.fq2} "
        "{input.fq1} {input.fq2} > {log}"

rule fastqc_clean_reads:
    input:
        fq1 = "output/2.ProcessedReads/{sample}_R1_clean.fastq.gz",
        fq2 = "output/2.ProcessedReads/{sample}_R2_clean.fastq.gz",
        directory_check = "output/1.QCs",
        directory_check2 = "output/3.QCsRmAdaptors"
    output: "output/0.LOGs/logPostQCs_{sample}.log"
    threads: 1
    message: "QC of trimmed reads on sample {wildcards.sample} files {input.fq1} {input.fq2}"
    conda:
        "10xmethylomes"
    shell: "fastqc --quiet -o output/3.QCsRmAdaptors/ {input.fq1} {input.fq2} > {output}"


# ---- MAPPING ---- # 

rule bowtie2_map:
    input:
        fq1 = "output/2.ProcessedReads/{sample}_R1_clean.fastq.gz",
        fq2 = "output/2.ProcessedReads/{sample}_R2_clean.fastq.gz"
    output: "output/4.Alignment/mapped/{sample}.bam"
    message: "Mapping of sample {wildcards.sample} files {input.fq1} {input.fq2}"
    log: "output/0.LOGs/logBowtie_{sample}.log"
    conda:
        "bulkatac"
    threads: 14
    shell:
        "(bowtie2 -q --threads 12 -X2000 -x {index} -1 {input.fq1} -2 {input.fq2} | "
        "samtools view -bSu - | "
        "samtools sort -T {wildcards.sample}_sorted - > {output}) > {log}"

rule flagStats_mapped:
    input: "output/4.Alignment/mapped/{sample}.bam"
    output: "output/4.Alignment/mapped/{sample}.flagStats"
    message: "Flagstat mapped {wildcards.sample}"
    threads: 1
    conda:
        "bulkatac"
    shell: "samtools flagstat {input} > {output}"


# ---- BAM FILE PROCESSING ---- # 

rule properPairs:
    input: "output/4.Alignment/mapped/{sample}.bam"
    output: temp("output/4.Alignment/mapped_pairs/{sample}_pairs.bam")
    message: "Output proper pairs for {wildcards.sample}"
    conda:
        "bulkatac"
    threads: 12
    shell: "/home/sbuckberry/working_data_01/bin/sambamba_v0.6.3 view -t {threads} -f bam -F 'proper_pair' {input} -o {output}"

rule flagStats_pairs:
    input: "output/4.Alignment/mapped_pairs/{sample}_pairs.bam"
    output: "output/4.Alignment/mapped_pairs/{sample}_pairs.flagStats"
    message: "Flagstat pairs {wildcards.sample}"
    threads: 1
    conda:
        "bulkatac"
    shell: "samtools flagstat {input} > {output}"

rule deduplicate:
    input: "output/4.Alignment/mapped_pairs/{sample}_pairs.bam"
    output: temp("output/4.Alignment/mapped_pairs_dedup/{sample}_pairs_dedup.bam")
    message: "Deduplicate reads of {wildcards.sample}"
    log: "output/0.LOGs/logDedup_{sample}.log"
    conda:
        "bulkatac"
    threads: 12
    shell: "/home/sbuckberry/working_data_01/bin/sambamba_v0.6.3 markdup -p -t {threads} --hash-table-size=1000000 {input} {output} 2> {log}"

rule flagStats_dedup:
    input: "output/4.Alignment/mapped_pairs_dedup/{sample}_pairs_dedup.bam"
    output: "output/4.Alignment/mapped_pairs_dedup/{sample}_pairs_dedup.flagStats"
    message: "Flagstat dedup {wildcards.sample}"
    threads: 1
    conda:
        "bulkatac"
    shell: "samtools flagstat {input} > {output}"

rule filter_blacklistedRegions:
    input: "output/4.Alignment/mapped_pairs_dedup/{sample}_pairs_dedup.bam"
    output:
        proper = temp("output/4.Alignment/mapped_pairs_dedup_filt/{sample}_pairs_dedup_filt.bam"),
        crap = temp("output/4.Alignment/mapped_pairs_dedup_filt/{sample}_pairs_dedup_filt.crap")
    message: "Filter out blacklisted regions for {wildcards.sample}"
    threads: 1
    conda:
        "bulkatac"
    shell: "samtools view -q 30 -L {blacklisted} -U {output.proper} -b {input} > {output.crap}"

rule flagStats_filt:
    input: "output/4.Alignment/mapped_pairs_dedup_filt/{sample}_pairs_dedup_filt.bam"
    output: "output/4.Alignment/mapped_pairs_dedup_filt/{sample}_pairs_dedup_filt.flagStats"
    message: "Flagstat for filtered {wildcards.sample}"
    threads: 1
    shell: "samtools flagstat {input} > {output}"

rule indx_filt:
    input: "output/4.Alignment/mapped_pairs_dedup_filt/{sample}_pairs_dedup_filt.bam"
    output: temp("output/4.Alignment/mapped_pairs_dedup_filt/{sample}_pairs_dedup_filt.bam.bai")
    message: "Index the intermitent (filtered) bam file for {wildcards.sample}"
    threads: 1
    conda:
        "bulkatac"
    shell: "samtools index {input}"

rule count_mitoReads:
    input: 
        bam = "output/4.Alignment/mapped_pairs_dedup_filt/{sample}_pairs_dedup_filt.bam",
        bai = "output/4.Alignment/mapped_pairs_dedup_filt/{sample}_pairs_dedup_filt.bam.bai"
    output: "output/4.Alignment/mapped_pairs_dedup_filt_mitocounts/{sample}_pairs_dedup_filt.mitoReads"
    message: "Count the mito reads of {wildcards.sample}"
    threads: 1
    conda:
        "bulkatac"
    shell: "samtools view -c {input.bam} chrM > {output}"

rule filter_mitoReads:
    input: "output/4.Alignment/mapped_pairs_dedup_filt/{sample}_pairs_dedup_filt.bam"
    output:
        proper = "output/4.Alignment/mapped_pairs_dedup_filt_noMT/{sample}_pairs_dedup_filt_noMT.bam",
        crap = temp("output/4.Alignment/mapped_pairs_dedup_filt_noMT/{sample}_pairs_dedup_filt_noMT.crap")
    message: "Remove the mito reads of {wildcards.sample}"
    threads: 1
    conda:
        "bulkatac"
    shell: "samtools view -q 30 -L {blacklistedMito} -U {output.proper} -b {input} > {output.crap}"

rule flagStats_noMT:
    input: "output/4.Alignment/mapped_pairs_dedup_filt_noMT/{sample}_pairs_dedup_filt_noMT.bam"
    output: "output/4.Alignment/mapped_pairs_dedup_filt_noMT/{sample}_pairs_dedup_filt_noMT.flagStats"
    message: "Flagstat of the final bam {wildcards.sample}"
    threads: 1
    conda:
        "bulkatac"
    shell: "samtools flagstat {input} > {output}"

rule indx_noMT:
    input: "output/4.Alignment/mapped_pairs_dedup_filt_noMT/{sample}_pairs_dedup_filt_noMT.bam"
    output: "output/4.Alignment/mapped_pairs_dedup_filt_noMT/{sample}_pairs_dedup_filt_noMT.bam.bai"
    message: "Index the final bam file for {wildcards.sample}"
    threads: 1
    conda:
        "bulkatac"
    shell: "samtools index {input}"


# ---- PEAKS AND COVERAGE ---- # 

rule coverage: 
    input:
        bam = "output/4.Alignment/mapped_pairs_dedup_filt_noMT/{sample}_pairs_dedup_filt_noMT.bam",
        bai = "output/4.Alignment/mapped_pairs_dedup_filt_noMT/{sample}_pairs_dedup_filt_noMT.bam.bai"
    output: "output/4.Alignment/coverage/{sample}_pairs_dedup_filt_noMT.coverage.bigwig"
    message: "Create bigwig files for {wildcards.sample}"
    conda:
        "deeptoolsenv"
    threads: 10
    shell: "bamCoverage --binSize 1 -p {threads} -b {input.bam} -o {output}"

rule readsInPromoters:
    input:
        bam = "output/4.Alignment/mapped_pairs_dedup_filt_noMT/{sample}_pairs_dedup_filt_noMT.bam",
        bai = "output/4.Alignment/mapped_pairs_dedup_filt_noMT/{sample}_pairs_dedup_filt_noMT.bam.bai"
    output: "output/4.Alignment/readsInProm/{sample}_pairs_dedup_filt_noMT.ReadsInProm"
    message: "Count the reads overlapping promoters for {wildcards.sample}"
    threads: 1
    conda:
        "bulkatac"
    shell: "samtools view -c -L {promoters} -b {input.bam} > {output}"

rule macs2PeakCallingNFR:
    input:
        bam = "output/4.Alignment/mapped_pairs_dedup_filt_noMT/{sample}_pairs_dedup_filt_noMT.bam",
        bai = "output/4.Alignment/mapped_pairs_dedup_filt_noMT/{sample}_pairs_dedup_filt_noMT.bam.bai"
    output: "output/4.Alignment/Peaks/{sample}_pairs_dedup_filt_noMT.NFR_peaks.narrowPeak"
    message: "Call NFR peaks for {wildcards.sample}"
    log: "output/0.LOGs/logMACS2NFR_{sample}.log"
    threads: 1
    params:
        nm = "output/4.Alignment/Peaks/{sample}_pairs_dedup_filt_noMT.NFR"
    conda:
        "deeptoolsenv"
    shell: "macs2 callpeak --nomodel --extsize 150 --shift -75 -t {input.bam} -f BAM -n {params.nm} --keep-dup all --gsize {size} 2> {log}"

rule macs2PeakCallingBAMPE: 
    input:
        bam = "output/4.Alignment/mapped_pairs_dedup_filt_noMT/{sample}_pairs_dedup_filt_noMT.bam",
        bai = "output/4.Alignment/mapped_pairs_dedup_filt_noMT/{sample}_pairs_dedup_filt_noMT.bam.bai"
    output: "output/4.Alignment/Peaks/{sample}_pairs_dedup_filt_noMT.BAMPE_peaks.narrowPeak"
    message: "Call BAMPE peaks for {wildcards.sample}"
    log: "output/0.LOGs/logMACS2BAMPE_{sample}.log"
    threads: 1
    params:
        nm = "output/4.Alignment/Peaks/{sample}_pairs_dedup_filt_noMT.BAMPE"
    conda:
        "deeptoolsenv"
    shell: "macs2 callpeak -t {input.bam} -f BAMPE -n {params.nm} --keep-dup all --gsize {size} 2> {log}"

rule NFR_readsInPeaks:
    input:
        bam = "output/4.Alignment/mapped_pairs_dedup_filt_noMT/{sample}_pairs_dedup_filt_noMT.bam",
        nfr_Peaks = "output/4.Alignment/Peaks/{sample}_pairs_dedup_filt_noMT.NFR_peaks.narrowPeak"
    output:
        rip = "output/4.Alignment/Peaks/{sample}_pairs_dedup_filt_noMT.ReadsInNFRpeaks",
        ripop = "output/4.Alignment/Peaks/{sample}_pairs_dedup_filt_noMT.ReadsInNFRpeaks_overProm"
    message: "Count reads overlapping NFR peaks for {wildcards.sample}"
    threads: 3
    conda:
        "bulkatac"
    shell:
        "cut -f1-3 {input.nfr_Peaks} | samtools view -L - -c -b {input.bam} > {output.rip};"
        "cut -f1-3 {input.nfr_Peaks} | /home/dmakosa/working_data_01/apps/miniconda3/envs/deeptoolsenv/bin/intersectBed -u -wa -a - -b {promoters} | samtools view -L - -c -b {input.bam} > {output.ripop}"

rule BAMPE_readsInPeaks:
    input:
        bam = "output/4.Alignment/mapped_pairs_dedup_filt_noMT/{sample}_pairs_dedup_filt_noMT.bam",
        bampe_Peaks = "output/4.Alignment/Peaks/{sample}_pairs_dedup_filt_noMT.BAMPE_peaks.narrowPeak"
    output:
        rip = "output/4.Alignment/Peaks/{sample}_pairs_dedup_filt_noMT.ReadsInBAMPEpeaks",
        ripop = "output/4.Alignment/Peaks/{sample}_pairs_dedup_filt_noMT.ReadsInBAMPEpeaks_overProm"
    message: "Count reads overlapping BAMPE peaks for {wildcards.sample}"
    threads: 3
    conda:
        "bulkatac"
    shell:
        "cut -f1-3 {input.bampe_Peaks} | samtools view -L - -c -b {input.bam} > {output.rip};"
        "cut -f1-3 {input.bampe_Peaks} | /home/dmakosa/working_data_01/apps/miniconda3/envs/deeptoolsenv/bin/intersectBed -u -wa -a - -b {promoters} | samtools view -L - -c -b {input.bam} > {output.ripop}"

rule generalPeakLandscape_readsInPeaks:
    input: "output/4.Alignment/mapped_pairs_dedup_filt_noMT/{sample}_pairs_dedup_filt_noMT.bam"
    output:
        rip = "output/4.Alignment/GeneralPeakLandscape/{sample}_pairs_dedup_filt_noMT.ReadsInGeneralPeakLandscape",
        ripop = "output/4.Alignment/GeneralPeakLandscape/{sample}_pairs_dedup_filt_noMT.ReadsInGeneralPeakLandscape_overProm"
    message: "Count reads overlapping peaks called on gold standard sample for {wildcards.sample}"
    threads: 3
    conda:
        "bulkatac"
    shell:
        "cut -f1-3 {generalPeakLanscape} | samtools view -L - -c -b {input} > {output.rip};"
        "cut -f1-3 {generalPeakLanscape} | /home/dmakosa/working_data_01/apps/miniconda3/envs/deeptoolsenv/bin/intersectBed -u -wa -a - -b {promoters} | samtools view -L - -c -b {input} > {output.ripop}"

rule OliviersPeakLanscape_readsInPeaks:
    input: "output/4.Alignment/mapped_pairs_dedup_filt_noMT/{sample}_pairs_dedup_filt_noMT.bam"
    output:
        rip = "output/4.Alignment/OliviersPeakLanscape/{sample}_pairs_dedup_filt_noMT.ReadsInOliviersPeakLanscape",
        ripop = "output/4.Alignment/OliviersPeakLanscape/{sample}_pairs_dedup_filt_noMT.ReadsInOliviersPeakLanscape_overProm"
    message: "Count reads overlapping peaks called on gold standard sample for {wildcards.sample}"
    threads: 3
    conda:
        "bulkatac"
    shell:
        "cut -f1-3 {OliviersPeakLanscape} | samtools view -L - -c -b {input} > {output.rip};"
        "cut -f1-3 {OliviersPeakLanscape} | /home/dmakosa/working_data_01/apps/miniconda3/envs/deeptoolsenv/bin/intersectBed -u -wa -a - -b {promoters} | samtools view -L - -c -b {input} > {output.ripop}"


# ---- CALCULATE FRIPs ---- # 

rule calculateNumberOfReads:
    input: "output/4.Alignment/mapped_pairs_dedup_filt_noMT/{sample}_pairs_dedup_filt_noMT.flagStats"
    output: "output/4.Alignment/mapped_pairs_dedup_filt_noMT/{sample}.noOfReads"
    message: "Calculate number of reads for {wildcards.sample}"
    threads: 2
    shell: "grep QC {input} | sed 's/ .*//g' > {output}"

shell.executable('/bin/bash')

rule calculateFRIPs:
    input:
        NoOfReadsOverPeaks = "output/4.Alignment/GeneralPeakLandscape/{sample}_pairs_dedup_filt_noMT.ReadsInGeneralPeakLandscape",
        NoOfReads = "output/4.Alignment/mapped_pairs_dedup_filt_noMT/{sample}.noOfReads"
    output: "output/4.Alignment/GeneralPeakLandscape/{sample}_pairs_dedup_filt_noMT.FRIP"
    message: "Calculate number of reads for {wildcards.sample}"
    threads: 6
    shell: "cat <(echo -e 'RL\tFRIP') "
            "<(paste <(echo {wildcards.sample}) "
            "<(paste <(cat {input.NoOfReadsOverPeaks}) <(echo '/') <(cat {input.NoOfReads}) <(echo '*100') | bc -l ) > {output})"

rule calculateFRIPs_OliviersLandscape:
    input:
        NoOfReadsOverPeaks = "output/4.Alignment/OliviersPeakLanscape/{sample}_pairs_dedup_filt_noMT.ReadsInOliviersPeakLanscape",
        NoOfReads = "output/4.Alignment/mapped_pairs_dedup_filt_noMT/{sample}.noOfReads"
    output: "output/4.Alignment/OliviersPeakLanscape/{sample}_pairs_dedup_filt_noMT.FRIPOlivier"
    message: "Calculate number of reads for {wildcards.sample}"
    threads: 6
    shell: "cat <(echo -e 'RL\tFRIP') "
            "<(paste <(echo {wildcards.sample}) "
            "<(paste <(cat {input.NoOfReadsOverPeaks}) <(echo '/') <(cat {input.NoOfReads}) <(echo '*100') | bc -l ) > {output})"

rule calculateFRIPromoters: # Add this one in to the file for re-computation
    input:
        NoOfReadsOverPromoters = "output/4.Alignment/readsInProm/{sample}_pairs_dedup_filt_noMT.ReadsInProm",
        NoOfReads = "output/4.Alignment/mapped_pairs_dedup_filt_noMT/{sample}.noOfReads"
    output: "output/4.Alignment/readsInProm/{sample}_pairs_dedup_filt_noMT.FRIPromoters"
    message: "Calculate percentage of reads in promoters for {wildcards.sample}"
    threads: 6
    shell: "cat <(echo -e 'RL\tFRIP') "
            "<(paste <(echo {wildcards.sample}) "
            "<(paste <(cat {input.NoOfReadsOverPromoters}) <(echo '/') <(cat {input.NoOfReads}) <(echo '*100') | bc -l ) > {output})"

# ---- CONCATENATE BAM FLAGSTATs ---- # 

rule combineBamFlagstats:
    input: 
        mapped = "output/4.Alignment/mapped/{sample}.flagStats",
        pairs = "output/4.Alignment/mapped_pairs/{sample}_pairs.flagStats",
        dedup = "output/4.Alignment/mapped_pairs_dedup/{sample}_pairs_dedup.flagStats",
        filt = "output/4.Alignment/mapped_pairs_dedup_filt/{sample}_pairs_dedup_filt.flagStats",
        nomt = "output/4.Alignment/mapped_pairs_dedup_filt_noMT/{sample}_pairs_dedup_filt_noMT.flagStats"
    output: "output/4.Alignment/flagstatsCombined/{sample}.flagStats_combined"
    message: "Combine FlagStats per sample for {wildcards.sample}"
    threads: 6
    shell: "cat <(echo -e 'RL\tTotalReads\tMapped\tPaired\tDedup\tFilt\tNoMT') "
            "<(paste <(echo {wildcards.sample}) "
            "<(grep QC {input.mapped} | sed 's/ .*//g') "
            "<(grep 'mapped (' {input.mapped} | sed 's/ .*//g') "
            "<(grep QC {input.pairs} | sed 's/ .*//g') "
            "<(grep QC {input.dedup} | sed 's/ .*//g') "
            "<(grep QC {input.filt} | sed 's/ .*//g') "
            "<(grep QC {input.nomt} | sed 's/ .*//g') ) > {output}"


# ---- INSERT SIZE PLOTS ---- # 

rule insertSize:
    input: "output/4.Alignment/mapped_pairs_dedup_filt_noMT/{sample}_pairs_dedup_filt_noMT.bam"
    output: "output/4.Alignment/InsertSizes/{sample}_pairs_dedup_filt_noMT.InsertSizesBAMPE"
    message: "Estimate insert size for {wildcards.sample}"
    threads: 3
    conda:
        "bulkatac"
    shell: "samtools view -f66 {input} | cut -f 9 | sed 's/^-//' > {output}"

rule insertSizePlot:
    input: 
        "output/4.Alignment/InsertSizes/{sample}_pairs_dedup_filt_noMT.InsertSizesBAMPE"
    output: "output/4.Alignment/InsertSizes/{sample}_pairs_dedup_filt_noMT.InsertSizesBAMPE_hist.png"
    threads: 1
    message: "Plot distribution of insert sizes for {wildcards.sample}"
    shell: "Rscript {histogram_script} {input} {output}"


# ---- SUBSAMPLING OF FILTERED BAM FILE ---- # 

rule calculateProportions:
    input: 
        numberOfReads = "output/4.Alignment/mapped_pairs_dedup_filt_noMT/{sample}.noOfReads",
        script = "/home/dmakosa/working_data_04/Bulk_ATAC_QCmetrics_tool/snakemake/helper_calculateProportions.sh"
    output: "output/5.SubsamplingOfFiltered/0.proportions/{sample}.proportions"
    message: "Calculate subsampling proportions for {wildcards.sample}"
    threads: 4
    shell: "sh {input.script} {input.numberOfReads} {wildcards.sample} {output}"

rule subsampleFilteredBAM:
    input:
        script = "/home/dmakosa/working_data_04/Bulk_ATAC_QCmetrics_tool/snakemake/helper_subsampleFilteredBAM.sh",
        proportions = "output/5.SubsamplingOfFiltered/0.proportions/{sample}.proportions"
    output: "output/5.SubsamplingOfFiltered/1.subsampling/subsampled_{sample}.txt"
    message: "Subsample filtered bam file for {wildcards.sample}"
    threads: 4 
    conda: "bulkatac"
    shell: "sh {input.script} {input.proportions} {wildcards.sample} {output} {threads}"

rule subsampleCallPeaks:
    input:
        script = "/home/dmakosa/working_data_04/Bulk_ATAC_QCmetrics_tool/snakemake/helper_subsampleCallPeaks.sh",
        subsamplingComplete = "output/5.SubsamplingOfFiltered/1.subsampling/subsampled_{sample}.txt",
        proportions = "output/5.SubsamplingOfFiltered/0.proportions/{sample}.proportions"
    output: "output/5.SubsamplingOfFiltered/2.callpeaks/NFRPeaks_{sample}.txt"
    message: "Subsample filtered bam file for {wildcards.sample}"
    threads: 2
    conda: "deeptoolsenv"
    shell: "sh {input.script} {input.proportions} {wildcards.sample} {output} {size}"

rule subsampleReadsInPeaks:
    input:
        script = "/home/dmakosa/working_data_04/Bulk_ATAC_QCmetrics_tool/snakemake/helper_subsampleReadsInPeaks.sh",
        proportions = "output/5.SubsamplingOfFiltered/0.proportions/{sample}.proportions",
        callingPeaksComplete = "output/5.SubsamplingOfFiltered/2.callpeaks/NFRPeaks_{sample}.txt"
    output: "output/5.SubsamplingOfFiltered/3.readsInPeaks/ReadsInPeaks_{sample}.txt"
    message: "Count reads overlapping NFR peaks for subsampled {wildcards.sample}"
    threads: 4
    conda: "bulkatac"
    shell: "sh {input.script} {input.proportions} {wildcards.sample} {output} {promoters}"

rule subsampleNumberOfReads:
    input:
        script = "/home/dmakosa/working_data_04/Bulk_ATAC_QCmetrics_tool/snakemake/helper_subsampleNumberOfReads.sh",
        proportions = "output/5.SubsamplingOfFiltered/0.proportions/{sample}.proportions",
        subsamplingComplete = "output/5.SubsamplingOfFiltered/1.subsampling/subsampled_{sample}.txt"
    output: "output/5.SubsamplingOfFiltered/3.readsInPeaks/NumberOfReads_{sample}.txt"
    threads: 2
    message: "Calculate number of reads for subsampled {wildcards.sample}"
    shell: "sh {input.script} {input.proportions} {wildcards.sample} {output}"

rule subsampleSummary:
    input:
        expand("output/5.SubsamplingOfFiltered/3.readsInPeaks/NumberOfReads_{sample}.txt", sample=config["samples"]),
        expand("output/5.SubsamplingOfFiltered/3.readsInPeaks/ReadsInPeaks_{sample}.txt", sample=config["samples"])
    output: "output/5.SubsamplingOfFiltered/4.summary/summary.tsv"
    message: "Summarize subsampled samples"
    threads: 6
    shell: "cat <(echo -e 'SampleNOR\tNoOfReads\tSampleRIP\tReadsInPeaks') "
            "<(paste <(printf '%s\n' output/5.SubsamplingOfFiltered/3.readsInPeaks/sub*NoOfReads) "
            "<(cat output/5.SubsamplingOfFiltered/3.readsInPeaks/sub*NoOfReads) "
            "<(printf '%s\n' output/5.SubsamplingOfFiltered/3.readsInPeaks/sub*ReadsInNFRpeaks) "
            "<(cat output/5.SubsamplingOfFiltered/3.readsInPeaks/sub*ReadsInNFRpeaks)) > {output}"

