#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : sourmash.snakefile
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 08.04.2020


rule solidify_sample:
    input:
        lambda wildcards: get_sample_files(wildcards.sample_accession)
    output:
        pipe('results/{sample_accession}/solid.{consumer}.pipe.fastq')
    resources:
        mem = 4000,
        time = lambda _: as_minutes(hours=16)
    threads: 1
    log:
        'logs/{sample_accession}/solidify.{consumer}.log'
    shell:
        'mkdir -p results/{wildcards.sample_accession}/ && '
        'goetia solid-filter -i {input} --pairing-mode split --solid-threshold {config[solid_threshold]} '
        '-K {config[K]} -x 2e9 -N 4 -o {output} > {log} 2>&1'


rule solid_sourmash:
    input:
        'results/{sample_accession}/solid.sourmash-stream.pipe.fastq'
    output:
        'results/{sample_accession}/solid.sourmash-stream.sig.gz'
    resources:
        mem = 8000,
        time = lambda _: as_minutes(hours=36)
    log:
        'logs/{sample_accession}/solid.sourmash-stream.log'
    shell:
        'goetia sourmash -N {config[minhash_size]} -K {config[K]} --pairing-mode single '
        '--names {wildcards.sample_accession}.solid -i {input} --save-sig /dev/stdout --save-stream '
        '--echo /dev/stderr 2> {log} | gzip --rsyncable -c > {output}'


rule solid_scaled_sourmash:
    input:
        'results/{sample_accession}/solid.sourmash-scaled-stream.pipe.fastq'
    output:
        'results/{sample_accession}/solid.sourmash-scaled-stream.sig.gz'
    resources:
        mem = 24000,
        time = lambda _: as_minutes(hours=36)
    log:
        'logs/{sample_accession}/solid.sourmash-scaled-stream.log'
    shell:
        'goetia sourmash --scaled {config[minhash_scaled]} -K {config[K]} --pairing-mode single '
        '--fine-interval 25000 '
        '--names {wildcards.sample_accession}.solid -i {input} --save-sig /dev/stdout --save-stream '
        '--echo /dev/stderr 2> {log} | gzip --rsyncable -c > {output}'


rule distances:
    input:
        'results/{sample_accession}/solid.{consumer}.sig.gz'
    output:
        'results/{sample_accession}/solid.{consumer}.distances.csv'
    resources:
        mem = lambda wildcards, attempt: attempt * 3000,
        time = lambda wildcards, attempt: as_minutes(hours=2**(attempt - 1))
    shell:
        'scripts/find_signature_distances.py {input} --output {output}'


rule solid_sourmash_all:
    input:
        expand('results/{sample_accession}/solid.{consumer}.distances.csv',
                sample_accession=samples.index.unique(),
                consumer=['sourmash-stream', 'sourmash-scaled-stream'])
