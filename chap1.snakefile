#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2022
# File   : common.snakefile
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 01.02.2022


GENOMIC_SAMPLES = read_sample_csv('config/chap1-genomes.csv')
GENOMIC_URLS = build_sample_urls(GENOMIC_SAMPLES)

TXOMIC_SAMPLES = read_sample_csv('config/chap1-txomes.csv')
TXOMIC_URLS = build_sample_urls(TXOMIC_SAMPLES)


rule download_stream_cdbg_build:
    conda: 'envs/goetia.yml'
    input:
        left  = 'data/streams/{accession}.pipe.1.fq.gz',
        right = 'data/streams/{accession}.pipe.2.fq.gz'
    output:
        'results/chap1/cdbg-stream/{accession}/goetia.cdbg.stats.json'
    threads: 3
    resources:
        mem = 32000,
        time = lambda _: as_minutes(hours=8)
    log: 'logs/cdbg-build/stream-{accession}.log'
    params:
        storage_type     = 'PHMapStorage',
        interval         = 5000000,
        K                = 31,
        track_cdbg_stats = True
    shell: '''
        goetia cdbg build -K {params.K} -S PHMapStorage -H FwdLemireShifter --interval {params.interval} \
        --track-cdbg-metrics --results-dir results/chap1/cdbg-stream/{wildcards.accession}/ --pairing-mode split \
        -i {input.left} {input.right}
    '''


rule results_figure_one:
    input:
        expand('results/chap1/cdbg-stream/{accession}/goetia.cdbg.stats.json',
               accession = GENOMIC_SAMPLES.index.unique())
