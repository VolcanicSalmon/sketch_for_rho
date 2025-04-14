import pandas as pd 
import glob
import pathlib
import os 
fws=glob.glob("*_R1.fastq.gz")
samp=[os.path.basename(file).split("_R1.fastq.gz")[0] for file in fws ]
rvs=[os.path.basename(file).replace("R1.fastq.gz","R2.fastq.gz") for file in fws]
rule all:
  input:
    "cm.txt"
rule rqc:
  input:
    fw="{samp}_R1.fastq.gz",
    rv="{samp}_R2.fastq.gz"
  output:
    fwhtml="{samp}_R1.fastqc.html",
    fwzip="{samp}_R1.fastqc.zip",
    rvhtml="{samp}_R2.fastqc.html",
    rvzip="{samp}_R2.fastqc.zip"
  threads: {}
  resources: {}
  shell:
    '''
    fastqc -t {threads} {input} 
    '''
rule trim:
  input: 
    fw="{samp}_R1.fastq.gz",
    rv="{samp}_R2.fastq.gz"
  output:
    fwtrim="{samp}_R1.trimmed.fq.gz",
    fwunp="{samp}_R1.unpaired.fq.gz",
    rvtrim="{samp}_R2.trimmed.fq.gz",
    rvunp="{samp}_R2.unpaired.fq.gz"
  params:
    minlen=config["trim_minlen"],
    minqual=config["trim_minqual"],
    adapter=config["trim_adapter"],
    winsize=config["trim_winsize"],
    threshold=config["trim_threshold"]
  shell:
    '''
    trimmomatic PE {input.fw} {input.rv} {output.fwtrim} {output.fwunp} {output.rvtrim} {output.rvunp} ILLUMINACLIP:{params.adapter} SLIDINGWINDOW:{params.winsize}:{params.threshold} LEADING:{params.minqual} TRAILING:{params.minqual} MINLEN:{params.minlen}
    '''
rule bwa_sam:
  input:
    fwtrim="{samp}_R1.trimmed.fq.gz",
    rvtrim="{samp}_R2.trimmed.fq.gz",
    ref=config["ref"]
  output:
    bam="{samp}.bam"
  params:
    mapthreads=contig["bwa_threads"]
  shell:
    '''
    bwa index {input.ref}
    bwa mem -t {params.mapthreads} {input.fwtrim} {input.rvtrim} |samtools view -hb -@ {params.mapthreads}| samtools sort -@ {params.mapthreads} - > {output.bam}
    '''
rule bam_depth:
  input:
    inbam="{samp}.bam",
    ref=config["ref"]
  output:
    coverage="{samp}.cover.txt"
  shell:
    '''
    samtools index {input.inbam}
    samtools faidx {input.ref} 
    samtools depth {input.inbam} | awk '{{all+=$3}} END {print "avg cov=",all/NR,NR}' > {output.coverage}
    '''
rule bam_dedup:
  input:
    inbam="{samp}.bam"
  output:
    outbam="{samp}_derep.bam"
    metrics="{samp}_metrics.txt"
  shell:
    '''
    gatk MarkDuplicates --INPUT {input.inbam} --METRICS_FILE {output.metrics} --OUTPUT "{output.outbam}"
    '''
rule bam_addrg:
  input:
    inbam="{samp}_derep.bam"
  output:
    outbam="{samp}_rg.bam",
    outbai="{samp}_rg.bai"
  params:
    sm="{samp}"
    pl=config["platform"]
  shell:
    '''
    gatk AddOrReplaceReadGroups -I {input.inbam} -O {output.outbam} --SORT_ORDER coordinate --RGSM {params.sm} --RGPU none --RGID {params.sm} --RGLB lib1 --RGPL lib --RGPL {params.pl}
    gatk BuildBamIndex I={output.outbam}
    '''
rule gatk_call:
  input:
    inbam="{samp}_rg.bam",
    ref=config["ref"]
  output:
    gvcf="{samp}_rg.g.vcf"
  params:
    mbq=config["bq"]
  shell:
    '''
    ls *dict 1> /dev/null 2>&1 ||   gatk CreateSequenceDictionary R={input.ref}
    gatk HaplotypeCaller -R {input.ref} -I {input.inbam} -O {output.gvcf} -ERC GVCF --mbq {params.mbq}
    '''
rule gatk_gendb:
  input:
    inbam="{samp}_rg.bam",
    ref=config["ref"],
    metadata=config["sample_metadata"]
  output:
    directory("mygendb/{interval}")
  params:
    interval="interval.txt"
    tmpdir=config["gatk_tempdb"]
    batchsize=config["gatk_batchsize"]
  shell:
    '''
    [ ! -f interval.txt ] && awk '{print $1}' {input.ref} > interval.txt
    gatk --java-options {params.java} -L {params.interval} --merge-input-intervals --tmp-dir {params.tmpdir} --batch-size {params.batchsize} --sample-name-map {input.metadata} --genomicsdb-workspace-path {output.gendb}
    '''
rule gatk_genotype:
  input:
    interval="interval.txt",
    ref=config["ref"],
  output:
    combinedvcf="combined.vcf.gz",
  params:
    gendb="mygendb"
    java=config["java"]
    vcf="{samp}.gt.vcf.gz"
  shell:
    '''
    gatk --java-options {params.java} GenotypeGVCFs -R {input.ref} -O {params.vcf} -V gendb://{params.gendb} -L {input.intervals} 
    vcfarray=$(ls *.gt.vcf.gz | awk '{{print "I="$1}}')
    gatk MergeVcfs $vcfarray -O {output.combinedvcf}
    '''
rule vcf_ped:
  input:
    invcf="combined.vcf.gz",
    ped=config["pedigree"]
  output:
    outvcf="combined_ped.vcf.gz"
  shell:
    '''
    gatk CalculateGenotypePosteriors -V {input.invcf} -ped {input.ped} -O {output.outvcf} 
    '''
rule lenient_ft:
  input:
    invcf="combined_ped.vcf.gz"
  output:
    vcf_rmis="combined_rmis.vcf.gz"
  params:
    fmis="F_MISSING > 0"

  shell:
    '''
    bcftools view -e {params.fmis} {input.invcf} -Oz -o {output.vcf_rmis} 
    '''
rule vcf_phase:
  input:
    invcf="combined_rmis.vcf.gz",
    ref=config["ref"],
    bam="{samp}_derep.bam"
  output:
    outphase="phased.vcf.gz"
  shell:
    '''
    whatshap phase --reference {input.ref} {input.invcf} {input.bam} -o {output.outphase}
    '''
rule vcf_ldhat:
  input:
    invcf="phased.vcf.gz"
  output:
    outlink="phased.link"
  params:
    lift_anno=config["lifted_chr"]

  shell:
    '''
    vcftools --gzvcf {input.invcf} --out {output.outlink} --ldhat --phased
    '''
rule vcf_likelihood:
  input:
    linkin="phased.link"
  output:
    liketab="phased.lik.tsv"
  params:
    ped=config["pedigree"]
  shell:
    '''
    vcf_to_lik.py -i {input.linkin} -o {output.liketab} -p {params.ped}
    '''
rule ldhat_interval:
  input:
    linkin=

