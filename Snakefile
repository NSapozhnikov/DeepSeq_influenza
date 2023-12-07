'''
Snakemake pipeline to process deep sequencing data
Steps:
-
-
-
'''

import os


sample_ids = ['SRR1705851', 'SRR1705858', 'SRR1705859', 'SRR1705860']

rule all:
    input:
        expand(os.path.join('data', 'SRR1705851_ref.fasta.amb'), sample=sample_ids),
        expand(os.path.join('data', 'SRR1705851_ref.fasta.ann'), sample=sample_ids),
        expand(os.path.join('data', 'SRR1705851_ref.fasta.bwt'), sample=sample_ids),
        expand(os.path.join('data', 'SRR1705851_ref.fasta.pac'), sample=sample_ids),
        expand(os.path.join('data', 'SRR1705851_ref.fasta.sa'), sample=sample_ids),
        expand(os.path.join('data', 'SRR1705851_ref.{sample}.varscan_results.vcf'), sample=sample_ids)


# rule get_roommate_data:
#     output: 
#         os.path.join('data', '{sample}_seqs.fastq.gz')
#     shell:
#         'wget -O {output} http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/001/SRR1705851/'


rule get_seqs:
    output:
        os.path.join('data', '{sample}_seqs.fastq.gz')
    params:
        url=lambda wildcards: f'http://ftp.sra.ebi.ac.uk/vol1/fastq/{wildcards.sample[:6]}/00{wildcards.sample[-1]}/{wildcards.sample}/{wildcards.sample}.fastq.gz'
    shell:
        'wget -O {output} {params.url}'

rule bwa_index:
    input:
        reference=os.path.join('data', 'SRR1705851_ref.fasta')
    output:
        os.path.join('data', '{reference}.fasta.amb'),
        os.path.join('data', '{reference}.fasta.ann'),
        os.path.join('data', '{reference}.fasta.bwt'),
        os.path.join('data', '{reference}.fasta.pac'),
        os.path.join('data', '{reference}.fasta.sa')
    shell:
        'bwa index {input.reference}'


rule bwa_align:
    input:
        os.path.join('data', '{reference}.fasta.amb'),
        os.path.join('data', '{reference}.fasta.ann'),
        os.path.join('data', '{reference}.fasta.bwt'),
        os.path.join('data', '{reference}.fasta.pac'),
        os.path.join('data', '{reference}.fasta.sa'),
        ref=os.path.join('data', '{reference}.fasta'),
        reads=os.path.join('data', '{sample}_seqs.fastq.gz')
    threads: 12
    log:
        'logs/bwa.{reference}.{sample}.log'
    output:
        temporary(os.path.join('data', '{reference}.{sample}.unsorted.bam'))
    shell:
        'bwa mem -t {threads} {input.ref} {input.reads} 2>{log} | samtools view -b > {output}'


rule bam_sort:
    input:
        rules.bwa_align.output     
    output:
        protected(os.path.join('data', '{reference}.{sample}.sorted.bam'))
    threads: 12
    shell:
        'samtools sort --threads {threads} {input} > {output}'

rule index_sorted_bam:
    input:
        rules.bam_sort.output
    output:
        protected(os.path.join('data', '{reference}.{sample}.sorted.indexed.bam.bai'))
    threads: 12
    shell:
        'samtools index {input} {output} -@ {threads}'

rule mpileup_seqs:
    input:
        ref=os.path.join('data', '{reference}.fasta'),
        bam=rules.bam_sort.output
    output:
        os.path.join('data', '{reference}.{sample}.sorted.indexed.mpileup')
    threads: 12
    shell:
        'samtools mpileup -d 0 -f {input.ref} {input.bam} >  {output}'

rule variant_calling:
    input: 
        os.path.join('data', '{reference}.{sample}.sorted.indexed.mpileup')
    output: 
        os.path.join('data', '{reference}.{sample}.varscan_results.vcf')
    shell:
        'java -jar /mnt/d/Bioinformatics/VarScan/VarScan.v2.3.9.jar mpileup2snp {input} --min-var-freq 0.001 --variants --output-vcf 1 > {output}'