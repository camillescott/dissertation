#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2022
# File   : common.snakefile
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 01.02.2022

GENOMIC_SAMPLES = read_sample_csv('config/chap2-genomes.csv')
GENOMIC_URLS = build_sample_urls(GENOMIC_SAMPLES)

TXOMIC_SAMPLES = read_sample_csv('config/chap2-txomes.csv')
TXOMIC_URLS = build_sample_urls(TXOMIC_SAMPLES)


rule stream_sourmash:
    conda: 'envs/goetia.yml'
    input:
        left  = 'data/fastx/{accession}.1.fq.gz',
        right = 'data/fastx/{accession}.2.fq.gz'
    output:
        stream = 'results/chap2/{accession}/sourmash-stream/sigs.gz',
        final = 'results/chap2/{accession}/sourmash-stream/final.sig'
    resources:
        mem      = 8000,
        time     = lambda _: as_minutes(hours=36)
    params:
        interval = 1000000,
        K        = 31,
        N        = 50000
    log:
        'logs/chap2/{accession}/sourmash-stream/log'
    shell:
        'goetia sketch sourmash -N {params.N} -K {params.K} '
        '--pairing-mode split --names {wildcards.accession} '
        '-i {input.left} {input.right} --save-stream /dev/stdout --save-stream-tick 5 '
        '--interval  {params.interval} '
        '--save-sig {output.final} '
        '2> {log} | gzip --rsyncable -c > {output.stream}'


rule stream_sourmash_distances:
    conda: 'envs/goetia.yml'
    input:
        'results/chap2/{accession}/sourmash-stream/sigs.gz'
    output:
        'results/chap2/{accession}/sourmash-stream/distances.csv'
    resources:
        mem = lambda wildcards, attempt: attempt * 3000,
        time = lambda wildcards, attempt: as_minutes(hours=2**(attempt - 1))
    shell:
        'scripts/find_signature_distances.py {input} --output {output} --signature-type SourmashSignature'


rule stream_sourmash_ref_distances:
    conda: 'envs/goetia.yml'
    input:
        'results/chap2/{accession}/sourmash-stream/sigs.gz'
    output:
        'results/chap2/{accession}/sourmash-stream/ref-distances.csv'
    resources:
        mem = lambda wildcards, attempt: attempt * 3000,
        time = lambda wildcards, attempt: as_minutes(hours=2**(attempt - 1))
    shell:
        'scripts/find_signature_distances_from_ref.py {input} --output {output} --signature-type SourmashSignature'


rule stream_draff:
    conda: 'envs/goetia.yml'
    input:
        left  = 'data/fastx/{accession}.1.fq.gz',
        right = 'data/fastx/{accession}.2.fq.gz'
    output:
        stream = 'results/chap2/{accession}/draff-stream/sigs.gz',
        final = 'results/chap2/{accession}/draff-stream/final.sig'
    resources:
        mem = 8000,
        time = lambda _: as_minutes(hours=36)
    params:
        interval = 1000000,
        W        = 31,
        K        = 9
    log:
        'logs/chap2/{accession}/draff-stream/log'
    shell:
        'goetia sketch draff -W {params.W} -K {params.K} '
        '-S BitStorage -x 4e5 --interval {params.interval} '
        '--save-sig {output.final} --save-stream /dev/stdout --save-stream-tick 5 '
        '--pairing-mode split -i {input.left} {input.right} '
        '2> {log} | gzip --rsyncable -c > {output.stream}'


rule stream_draff_distances:
    conda: 'envs/goetia.yml'
    input:
        'results/chap2/{accession}/draff-stream/sigs.gz'
    output:
        'results/chap2/{accession}/draff-stream/distances.csv'
    resources:
        mem = lambda wildcards, attempt: attempt * 3000,
        time = lambda wildcards, attempt: as_minutes(hours=2**(attempt - 1))
    shell:
        'scripts/find_signature_distances.py {input} --output {output} '
        '--signature-type DraffSignature --distance-metric cosine canberra '
        'correlation braycurtis'


rule stream_draff_ref_distances:
    conda: 'envs/goetia.yml'
    input:
        'results/chap2/{accession}/draff-stream/sigs.gz'
    output:
        'results/chap2/{accession}/draff-stream/ref-distances.csv'
    resources:
        mem = lambda wildcards, attempt: attempt * 3000,
        time = lambda wildcards, attempt: as_minutes(hours=2**(attempt - 1))
    shell:
        'scripts/find_signature_distances_from_ref.py {input} --output {output} '
        '--signature-type DraffSignature --distance-metric cosine canberra '
        'correlation braycurtis'


rule all_sourmash_stream:
    input:
        expand('results/chap2/{accession}/sourmash-stream/final.sig',
               accession = list(TXOMIC_SAMPLES.index.unique()) + list(GENOMIC_SAMPLES.index.unique()))


rule all_sourmash_distances:
    input:
        expand('results/chap2/{accession}/sourmash-stream/distances.csv',
               accession = list(TXOMIC_SAMPLES.index.unique()) + list(GENOMIC_SAMPLES.index.unique()))


rule all_sourmash_ref_distances:
    input:
        expand('results/chap2/{accession}/sourmash-stream/ref-distances.csv',
               accession = list(TXOMIC_SAMPLES.index.unique()) + list(GENOMIC_SAMPLES.index.unique()))


rule all_draff_stream:
    input:
        expand('results/chap2/{accession}/draff-stream/final.sig',
               accession = list(TXOMIC_SAMPLES.index.unique()) + list(GENOMIC_SAMPLES.index.unique()))


rule all_draff_distances:
    input:
        expand('results/chap2/{accession}/draff-stream/distances.csv',
               accession = list(TXOMIC_SAMPLES.index.unique()) + list(GENOMIC_SAMPLES.index.unique()))


rule all_draff_ref_distances:
    input:
        expand('results/chap2/{accession}/draff-stream/ref-distances.csv',
               accession = list(TXOMIC_SAMPLES.index.unique()) + list(GENOMIC_SAMPLES.index.unique()))


rule all_distances:
    input:
        rules.all_sourmash_distances.input,
        rules.all_sourmash_ref_distances.input,
        rules.all_draff_distances.input,
        rules.all_draff_ref_distances.input,
