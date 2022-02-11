#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : data.snakefile
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 08.04.2020


rule download_runs:
    output:
        'data/runs/{filename}'
    params:
        url = lambda wildcards: samples[samples['fastq_filename'] == f'{wildcards.filename}'].fastq_ftp[0]
    resources:
        mem = 1000,
        time = lambda _: as_minutes(hours=16)
    log: 'logs/download/{filename}.download.log'
    message:
        'Download given run fastq. If --ignore-incomplete and --keep-incomplete are used, resumes unfinished download.'
    shell:
        'curl --connect-timeout 10 --keepalive-time 2 -C - -o {output} -L {params.url} > {log} 2>&1'


rule download:
    input:
        [f'data/runs/{filename}' for filename in samples.fastq_filename]


for sample in samples.index.unique():
    files = get_sample_files(sample)
    files = '\n\t- '.join(files)
    print(f'{sample}:\n\t- {files}')

