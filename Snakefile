import pandas as pd
import os

##GLOBALS##

configfile: "config.yaml"

raw=config["raw"]

samples = pd.read_table(config["samples"])
attributes=config["attributes"]
compare = pd.read_table(config["compare"])

reps = ["1","2"]

host_genome=config["host_genome"]
spike_genome=config["spike_genome"]
host_ucsc_len_file=config["host_ucsc_len_file"]
chrmap=config["chrmap"]

filters = ["total","unique","primary","unique_"+spike_genome,"unique_"+host_genome,"primary_"+spike_genome,"primary_"+host_genome]
report_filters = ["total","primary_"+host_genome,"unique_"+host_genome]
output_filters = ["unique_"+host_genome,"primary_"+host_genome]

bwa_index=config["bwa_index"]
adapter1=config["adapter1"]
adapter2=config["adapter2"]

##Functions

##Rules

rule all:
  input:
    expand("fastq/{sample}.{read}.fq.gz",sample=samples.Sample,read=[1,2]),
    expand("fastq/{sample}.{read}_fastqc.html",sample=samples.Sample,read=[1,2]),
    "fastq/multiqc_report.html",
    expand("cutadapt/{sample}.{read}.cutadapt.fq.gz",sample=samples.Sample,read=[1,2]),
    expand("cutadapt/{sample}.{read}.cutadapt_fastqc.html",sample=samples.Sample,read=[1,2]),
    expand("cutadapt/{sample}.{read}.cutadapt_screen.txt",sample=samples.Sample,read=[1,2]),
    "cutadapt/multiqc_report.html",
    expand("bwa/{sample}/{sample}.{filter}.flagstat",sample=samples.Sample,filter=report_filters),
    expand("bwa/{sample}/{sample}.{filter}.bam",sample=samples.Sample,filter=report_filters),
    expand("bwa/{sample}/{sample}.{filter}.bam.bai",sample=samples.Sample,filter=report_filters),
    "bwa/multiqc_report.html",
    "bwa/readStats.html",
    expand("visualisation_log2/{out_filters}/{IP}_{rep}.{out_filters}.input.log2.bpm.bw",IP=compare.IP,rep=reps,out_filters=output_filters),
    expand("visualisation_log2/{out_filters}/{IP}_{rep}.{out_filters}.input.log2.bpm.ucsc.bw",IP=compare.IP,rep=reps,out_filters=output_filters),
    expand("visualisation/{sample}.{out_filters}.bpm.bw",sample=samples.Sample,out_filters=output_filters),
    #"deeptools/unique_matrix.npz",
    #"deeptools/unique_correlation.pdf",
    #"deeptools/unique_pca.pdf",
    expand("macs2_unique/{IP}_{rep}_peaks.narrowPeak",IP=compare.IP,rep=reps),
    expand("macs2_primary/{IP}_{rep}_peaks.narrowPeak",IP=compare.IP,rep=reps),
    expand("macs2_unique/{IP}_{rep}_peaks.broadPeak",IP=compare.IP,rep=reps),
    expand("macs2_primary/{IP}_{rep}_peaks.broadPeak",IP=compare.IP,rep=reps)

rule get_fastq:
  input:
     raw+"/{sample}.{read}.fq.gz"
  output:
     "fastq/{sample}.{read}.fq.gz"
  shell:
    "ln -s {input} {output}"

rule fastqc_raw:
  input:
   "fastq/{sample}.{read}.fq.gz"
  threads: 4
  output:
   "fastq/{sample}.{read}_fastqc.html"
  shell:
   "fastqc-0.11.9 {input}"

rule multiqc_raw:
  input:
    expand("fastq/{sample}.{read}_fastqc.html",sample=samples.Sample,read=[1,2])
  output:
    "fastq/multiqc_report.html"
  shell:
    "multiqc -f fastq -o fastq"

rule cutadapt:
  input:
    r1="fastq/{sample}.1.fq.gz",
    r2="fastq/{sample}.2.fq.gz"
  params:
    a1=adapter1,
    a2=adapter2
  threads: 4
  output:
    p1="cutadapt/{sample}.1.cutadapt.fq.gz",
    p2="cutadapt/{sample}.2.cutadapt.fq.gz"
  shell:
    "cutadapt-1.18 -a {params.a1} -A {params.a2} -j {threads} --minimum-length 20 --nextseq-trim=20 -o {output.p1} -p {output.p2} {input.r1} {input.r2}"

rule cutadapt_fastqc:
  input:
   "cutadapt/{sample}.{read}.cutadapt.fq.gz"
  output:
   "cutadapt/{sample}.{read}.cutadapt_fastqc.html"
  shell:
   "fastqc-0.11.9 {input}"

rule cutadapt_fastq_screen:
  input:
    "cutadapt/{sample}.{read}.cutadapt.fq.gz"
  output:
    "cutadapt/{sample}.{read}.cutadapt_screen.txt"
  shell:
    "fastq_screen --conf /homes/genomes/tool_configs/fastq_screen/fastq_screen.conf --outdir cutadapt {input}"

rule cutadapt_multiqc:
  input:
    expand("cutadapt/{sample}.{read}.cutadapt_fastqc.html",sample=samples.Sample,read=[1,2]),
    expand("cutadapt/{sample}.{read}.cutadapt_screen.txt",sample=samples.Sample,read=[1,2])
  output:
    "cutadapt/multiqc_report.html"
  shell:
    "multiqc -f cutadapt -o cutadapt"

rule mapBWA:
  input:
    r1="cutadapt/{sample}.1.cutadapt.fq.gz",
    r2="cutadapt/{sample}.2.cutadapt.fq.gz"
  params:
    index=bwa_index
  output:
    temp("bwa/{sample}/{sample}.sam")
  threads: 5
  shell:
    "bwa mem -t {threads} -M {params.index} {input.r1} {input.r2}  > {output}"

rule processBam:
  input:
   "bwa/{sample}/{sample}.sam"
  threads: 5
  params:
   prefix="{sample}"
  output:
   "bwa/{sample}/{sample}.total.flagstat",
   "bwa/{sample}/{sample}.total.bam",
   "bwa/{sample}/{sample}.total.bam.bai",
   "bwa/{sample}/{sample}.primary.flagstat",
   "bwa/{sample}/{sample}.primary.bam",
   "bwa/{sample}/{sample}.primary.bam.bai",
   "bwa/{sample}/{sample}.unique.flagstat",
   "bwa/{sample}/{sample}.unique.bam",
   "bwa/{sample}/{sample}.unique.bam.bai",
   "bwa/{sample}/{sample}.primary_"+host_genome+".flagstat",
   "bwa/{sample}/{sample}.primary_"+host_genome+".bam",
   "bwa/{sample}/{sample}.primary_"+host_genome+".bam.bai",
   "bwa/{sample}/{sample}.unique_"+host_genome+".flagstat",
   "bwa/{sample}/{sample}.unique_"+host_genome+".bam",
   "bwa/{sample}/{sample}.unique_"+host_genome+".bam.bai"
  shell:
   "bash bin/scripts/processBam_paired.sh {params.prefix}"

rule readStats:
  input:
    expand("bwa/{sample}/{sample}.unique_"+host_genome+".flagstat",sample=samples.Sample)
  params:
    filt=",".join(report_filters),
    attr=attributes
  output:
    "bwa/readStats.html"
  shell:
    """
    R -e "rmarkdown::render('bin/scripts/readStats.Rmd')" --args "../../bwa" "{params.filt}" "../../samples.tsv" "{params.attr}"
    mv bin/scripts/readStats.html {output}
    """

rule bwa_multiqc:
  input:
    expand("bwa/{sample}/{sample}.{report_filt}.bam",sample=samples.Sample,report_filt=report_filters)
  params:
    files="unique_"+host_genome
  output:
    "bwa/multiqc_report.html"
  shell:
    "multiqc -f bwa/*/*{params.files}* -o bwa"

rule bamCoverage:
  input:
    "bwa/{sample}/{sample}.{out_filters}.bam"
  output:
    "visualisation/{sample}.{out_filters}.bpm.bw"
  threads: 12
  shell:
    """
    bamCoverage -bs 1 --normalizeUsing BPM -p {threads} -b {input} --outFileName {output}
    """

rule bamCompare:
  input:
    b1="bwa/{IP}_IP_{rep}/{IP}_IP_{rep}.{out_filters}.bam",
    b2="bwa/{IP}_Input_{rep}/{IP}_Input_{rep}.{out_filters}.bam"
  output:
    "visualisation_log2/{out_filters}/{IP}_{rep}.{out_filters}.input.log2.bpm.bw"
  threads: 4
  shell:
    """
    bamCompare -b1 {input.b1} -b2 {input.b2} -o {output} -bs 20 -p {threads} --normalizeUsing BPM --exactScaling --scaleFactorsMethod None
    """

rule bw_ucsc:
  input:
    "visualisation_log2/{out_filters}/{IP}_{rep}.{out_filters}.input.log2.bpm.bw"
  params:
    len=host_ucsc_len_file,
    chrmap=chrmap,
    bw="visualisation_log2/{out_filters}/{IP}_{rep}.{out_filters}.input.log2.bpm"
  output:
    "visualisation_log2/{out_filters}/{IP}_{rep}.{out_filters}.input.log2.bpm.ucsc.bw"
  shell:
    """
    bash bin/scripts/bw_ucsc.sh {params.bw} {params.chrmap} {params.len}
    """

rule macs2_peaks_unique:
  input:
    t="bwa/{IP}_IP_{rep}/{IP}_IP_{rep}.unique_"+host_genome+".bam",
    c="bwa/{IP}_Input_{rep}/{IP}_Input_{rep}.unique_"+host_genome+".bam"
  params:
    dir="macs2_unique/",
    out="{IP}_{rep}"
  output:
    "macs2_unique/{IP}_{rep}_peaks.narrowPeak",
    "macs2_unique/{IP}_{rep}_peaks.broadPeak"
  shell:
    """
    macs2 callpeak -t {input.t} -c {input.c} -g 1.1e9 -f BAMPE --outdir {params.dir} --name {params.out} --tempdir macs2_unique
    macs2 callpeak -t {input.t} -c {input.c} -g 1.1e9 -f BAMPE --outdir {params.dir} --name {params.out} --broad --tempdir macs2_unique
    """

rule macs2_peaks_primary:
  input:
    t="bwa/{IP}_IP_{rep}/{IP}_IP_{rep}.primary_"+host_genome+".bam",
    c="bwa/{IP}_Input_{rep}/{IP}_Input_{rep}.primary_"+host_genome+".bam"
  params:
    dir="macs2_primary/",
    out="{IP}_{rep}"
  output:
    "macs2_primary/{IP}_{rep}_peaks.narrowPeak",
    "macs2_primary/{IP}_{rep}_peaks.broadPeak"
  shell:
    """
    macs2 callpeak -t {input.t} -c {input.c} -g 1.1e9 -f BAMPE --outdir {params.dir} --name {params.out} --tempdir macs2_primary
    macs2 callpeak -t {input.t} -c {input.c} -g 1.1e9 -f BAMPE --outdir {params.dir} --name {params.out} --broad --tempdir macs2_primary
    """

### Need to rethink how to do these for each group of experiments?

rule multiBamSummary:
  input:
    expand("bwa/{sample}/{sample}.primary_"+host_genome+".bam",sample=samples.Sample)
  params:
    scale="deeptools/primary.sf.txt"
  threads: 20
  output:
    "deeptools/primary_matrix.npz"
  shell:
    """
    multiBamSummary-3.5.0 bins -b {input} -o {output} --smartLabels -p {threads} --scalingFactors {params.scale}
    """

rule plotCorrelation:
  input:
    "deeptools/primary_matrix.npz",
  output:
    uc="deeptools/primary_correlation.pdf",
    up="deeptools/primary_pca.pdf"
  shell:
    """
    plotCorrelation -in {input} -c pearson -p heatmap -o {output.uc} --removeOutliers
    plotPCA -in {input} -o {output.up}
    """

