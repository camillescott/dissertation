#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : data.snakefile
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 08.04.2020


ALL_SAMPLES = pd.concat([read_sample_csv('config/chap1-genomes.csv'),
                         read_sample_csv('config/chap1-txomes.csv'),
                         read_sample_csv('config/chap2-genomes.csv'),
                         read_sample_csv('config/chap2-txomes.csv')])
ALL_SAMPLES.drop_duplicates(inplace=True)
ACCESSIONS = ALL_SAMPLES.index.unique()

print('ACCESSIONS: ', ' '.join(ACCESSIONS))


rule download_sample_left:
    output:
        data = 'data/fastx/{accession}.1.fq.gz'
    params:
        url = lambda wildcards: paired_url_split(wildcards.accession,
                                                 ALL_SAMPLES)[0]
    shell: '''
        curl -sL -o {output.data} "{params.url}" 
    '''


rule download_sample_right:
    output:
        data = 'data/fastx/{accession}.2.fq.gz'
    params:
        url = lambda wildcards: paired_url_split(wildcards.accession,
                                                 ALL_SAMPLES)[1]
    shell: '''
        curl -sL -o {output.data} "{params.url}" 
    '''


rule download_stream_sample_left:
    output:
        data = pipe('data/streams/{accession}.pipe.1.fq.gz'),
        stats = 'results/chap1/{accession}.1.curl.txt'
    params:
        url = lambda wildcards: paired_url_split(wildcards.accession,
                                                 ALL_SAMPLES)[0]
    shell: '''
        set -o pipefail; curl -sL "{params.url}" | pv -tbn 2>> {output.stats} 1>> {output.data}
    '''


rule download_stream_sample_right:
    output:
        data = pipe('data/streams/{accession}.pipe.2.fq.gz'),
        stats = 'results/chap1/{accession}.2.curl.txt'
    params:
        url = lambda wildcards: paired_url_split(wildcards.accession,
                                                 ALL_SAMPLES)[1]
    shell: '''
        set -o pipefail; curl -sL "{params.url}" | pv -tbn 2>> {output.stats} 1>> {output.data}
    '''


rule download:
    input:
        expand('data/fastx/{accession}.{side}.fq.gz',
                side = [1,2],
                accession = ACCESSIONS)


#for sample in samples.index.unique():
#    files = get_sample_files(sample)
#    files = '\n\t- '.join(files)
#    print(f'{sample}:\n\t- {files}')

