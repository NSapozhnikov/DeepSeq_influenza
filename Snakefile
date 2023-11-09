'''
Snakemake pipeline to process deep sequencing data
Steps:
-
-
-
'''

import os


URL_seq = 'http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/001/SRR1705851/SRR1705851.fastq.gz'
# URL_ref = ''

rule all:
    input:
        os.path.join('data', 'SRR1705851_ref.SRR1705851_seqs.sorted.bam')


rule get_seqs:
    output: 
        os.path.join('data', 'SRR1705851_seqs.fastq.gz')
    shell:
        'wget -O {output} {URL_seq}'


# rule get_ref:
#     output:
#         os.path.join('data', 'SRR1705851_ref.fasta')
#     shell:
#         'wget -O {output} {URL_ref}'


rule bwa_index:
    input:
        os.path.join('data', '{reference}.fasta')
    output:
        os.path.join('data', '{reference}.fasta.amb'),
        os.path.join('data', '{reference}.fasta.ann'),
        os.path.join('data', '{reference}.fasta.bwt'),
        os.path.join('data', '{reference}.fasta.pac'),
        os.path.join('data', '{reference}.fasta.sa')
    shell:
        'bwa index {input}'


rule bwa_align:
    input:
        os.path.join('data', '{reference}.fasta.amb'),
        os.path.join('data', '{reference}.fasta.ann'),
        os.path.join('data', '{reference}.fasta.bwt'),
        os.path.join('data', '{reference}.fasta.pac'),
        os.path.join('data', '{reference}.fasta.sa'),
        ref=os.path.join('data', '{reference}.fasta'),
        reads=os.path.join('data', '{seqs}.fastq.gz')
    threads: 12
    log:
        'logs/bwa.{reference}.{seqs}.log'
    output:
        temporary(os.path.join('data', '{reference}.{seqs}.unsorted.bam'))
    shell:
        'bwa mem -t {threads} {input.ref} {input.reads} 2>{log} | samtools view -b > {output}'


rule bam_sort:
    input:
        rules.bwa_align.output     
    output:
        protected(os.path.join('data', '{reference}.{seqs}.sorted.bam'))
    threads: 12
    shell:
        'samtools sort --threads {threads} {input} > {output}'