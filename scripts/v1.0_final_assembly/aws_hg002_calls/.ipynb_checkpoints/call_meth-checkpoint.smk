#!/home/timp/miniconda3/bin/snakemake --snakefile
"""
This Snakemake pipeline is intended to absorb HG002 pangenome sequencing data and methylation call it on AWS.  Based on Ariel's original snakemake
    * Aligns to the reference genome with winnowmap
    * Indexes methylation with f5c
    * Calls methylation with Nanopolish
    * Formats methylation calls with nanopore-methylation-utilities
    fast5 directory and fastq file must have same basename, summary file not required 
"""
configfile: "config.yaml"
###-------input paths/files -------###

##WT - note check this to figure out parallelization properly
##cores=config["cores"]

##Software config
nanopolish = config["nanopolish"]
winnowmap = config["winnowmap"]
util = config["nanopore-methylation-utilities"]
f5c = config["f5c"]


##Pass in directory we are working in --config base=dir
base=config["base"]
##Pass in which references using --config ref=chm13 or ref=hg002 - default is chm13.
if config["ref"]=="hg002":
    ref = config["hg002-ref"]
else:
    ref = config["chm13-ref"]

###--------------###

##WT - again way more complicated than necessary. I noted that there is a pattern of GM24385_{0..22} I can use, broken only by {NB}

# ###------- Pipeline Rules -------#####
##Ok - this is supposed to basically request the final output, to make life simple
rule all:
    input:
        expand( base_out_path + "/meth_calls/{id}_meth.bam", sample1=SAMPLESLONG)
##Make reference index        
rule ref_index:
    output:
        ref+".k15.txt"
    threads: 96
    shell:
        """
        {winnowmap}/meryl count threads={threads} k=15 output merylDB {ref}
        {winnowmap}/meryl print greater-than distinct=0.9998 merylDB > {output}
        """
##Align with winnowmap
rule align:
    input:
        ##Question about how to "find" fastq
        ##Ok - it seems to me we don't *need* the fastq as an input paradoxically.  We just need the index as an input and we use the wildcard of the output to find the input?
        #fastq = fastqpath + "/{id}_Guppy_4.2.2_prom.fastq.gz",
        refidx=rules.ref_index.output
    output:
        bam = base + "/bam/{id}.bam"
    threads: 96
    message: """Aligning to reference with winnowmap"""
    run:
        ##Find fq block here
        shell("{winnowmap}/winnowmap -t {threads} -W {input.refidx} -ax map-ont {ref} {fq} | samtools view -b -u -F 256 | samtools sort -o {output.bam}")
        shell("samtools index {output.bam}")

#rule np_index:
##In principle I need this, but I already indexed so I'm not going to bother for now.
        
rule nanopolish:
    input:
        bam = rules.align.output.bam
    output:
        meth=base + "/meth_calls/{id}_CpG_methylation.tsv"
    params:
        DIR = base_out_path + "/meth_calls"
    message: """calling methylation with nanopolish"""
    threads: 96
    run:
        ##Find fq block here
    	shell("{nanopolish}/nanopolish call-methylation -b {input.bam} -r {fq} -g {ref} -q cpg -t {threads} --progress > {output.meth}")
    
rule methylbed:
    input:
        tsv = base + "/meth_calls/{id}_CpG_methylation.tsv"
    output:
        methbed=base + "/meth_calls/{sample1}_CpG_methylation.bed.gz"
    message: """format methyl bed"""
    threads: 1
    run:
        temp=base+"/meth_calls/{id}.tmp"
        shell("python3 {util}/mtsv2bedGraph.py -q cpg -c 1.5 -g {ref} -i {input.tsv} > {temp}")
        shell("sort {temp} -k1,1 -k2,2n | bgzip > {output.methbed}")
        shell("tabix -p bed {output.methbed}")
        shell("rm {temp}")

rule methylbam:
    input:
        bam = base + "/bam/{id}.bam",
        tsv = rules.methylbed.output.methbed
    output:
        base + "/meth_calls/{id}_meth.bam"
    params:
        DIR = base_out_path + "/meth_calls",
        file = "{sample1}_CpG_methylation.bed.gz",
        filt = base_out_path + "/bam/{sample1}_primary.bam",
        methbam = "{sample1}_meth.bam",
        ref = ref,
        cores = cores
    message: """convert bam for methylation"""
    threads: 96
    shell:
        """
        python3 {util}/convert_bam_for_methylation.py -t {threads} --verbose -b {input.bam} \
               -c {input.tsv} -f {ref} |\
                samtools sort -o {output}
        samtools index {output}
        """
        