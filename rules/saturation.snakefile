#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2021
# File   : sourmash.snakefile
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 22.10.2021


rule stream_diginorm:
    '''
    Digitally normalize sample runs to output pipe at specified coverage threshold.
    '''
    input:
        lambda wildcards: get_sample_files(wildcards.sample_accession)
    output:
        pipe('results/{sample_accession}/diginorm-{diginorm_coverage}/to-{consumer}.pipe.fastq')
    resources:
        mem = 8000,
        time = lambda _: as_minutes(hours=16)
    threads: 1
    log: 'logs/{sample_accession}/diginorm-{diginorm_coverage}/to-{consumer}.log'
    shell:
        'mkdir -p results/{wildcards.sample_accession}/diginorm-{diginorm_coverage} && ' 
        'goetia filter diginorm -i {input} --pairing-mode split '
        '-K {config[K]} -C {wildcards.diginorm_coverage} '
        '-S ByteStorage -H CanLemireShifter -x 1e9 -N 4 '
        '--interval 1000000 -o {output} > {log} 2>&1'


rule stream_merge:
    '''
    Merged paired reads into a single stream. This is technically a no-op, because
    the downstream programs can also read paired data, but it allows a no-filter
    producer without adding more downstream rules, and the parsing is ridiculously
    fast anyway.
    '''
    input:
        lambda wildcards: get_sample_files(wildcards.sample_accession)
    output:
        pipe('results/{sample_accession}/unfiltered/to-{consumer}.pipe.fastq')
    resources:
        mem = 8000,
        time = lambda _: as_minutes(hours=16)
    threads: 1
    log: 'logs/{sample_accession}/unfiltered/to-{consumer}.log'
    shell:
        'mkdir -p results/{wildcards.sample_accession}/unfiltered && ' 
        'goetia utils merge-paired -i {input} --pairing-mode split '
        '--interval 1000000 -o {output} > {log} 2>&1'


rule stream_solid:
    '''
    Filter out reads that don't have the specified proportion of k-mers
    with count over the given threshold.
    '''
    input:
        lambda wildcards: get_sample_files(wildcards.sample_accession)
    output:
        pipe('results/{sample_accession}/solid-{solid_min_count}/to-{consumer}.pipe.fastq')
    resources:
        mem = 8000,
        time = lambda _: as_minutes(hours=16)
    threads: 1
    log: 'logs/{sample_accession}/solid-{solid_min_count}/to-{consumer}.log'
    log: 'logs/{sample_accession}/diginorm-{diginorm_coverage}/to-{consumer}.log'
    shell:
        'mkdir -p results/{wildcards.sample_accession}/solid-{solid_min_count} && ' 
        'goetia filter solid -i {input} --pairing-mode split '
        '-K {config[K]} -C {wildcards.solid_min_count} -P 0.8 '
        '-S ByteStorage -H CanLemireShifter -x 1e9 -N 4 '
        '--interval 1000000 -o {output} > {log} 2>&1'


rule stream_sourmash:
    input:
        'results/{sample_accession}/{producer}/to-sourmash-stream.pipe.fastq'
    output:
        stream = 'results/{sample_accession}/sourmash-stream/from-{producer}.sigs.gz',
        final = 'results/{sample_accession}/sourmash-stream/from-{producer}.final.sig'
    resources:
        mem = 8000,
        time = lambda _: as_minutes(hours=36)
    log:
        'logs/{sample_accession}/sourmash-stream/from-{producer}.log'
    shell:
        'goetia sketch sourmash -N {config[minhash_size]} -K {config[K]} '
        '--pairing-mode single --names {wildcards.sample_accession}.{wildcards.producer} '
        '-i {input} --save-stream /dev/stdout --save-stream-tick 5 '
        '--interval {config[sketch_distance_interval]} '
        '--save-sig {output.final} '
        '2> {log} | gzip --rsyncable -c > {output.stream}'


rule stream_sourmash_distances:
    input:
        'results/{sample_accession}/sourmash-stream/from-{producer}.sigs.gz'
    output:
        'results/{sample_accession}/sourmash-stream/from-{producer}.stream-distances.csv'
    resources:
        mem = lambda wildcards, attempt: attempt * 3000,
        time = lambda wildcards, attempt: as_minutes(hours=2**(attempt - 1))
    shell:
        'scripts/find_signature_distances.py {input} --output {output} --signature-type SourmashSignature'



rule stream_sourmash_ref_distances:
    input:
        'results/{sample_accession}/sourmash-stream/from-{producer}.sigs.gz'
        #ref = 'results/{sample_accession}/sourmash-stream/from-{producer}.final.sig'
    output:
        'results/{sample_accession}/sourmash-stream/from-{producer}.ref-distances.csv'
    resources:
        mem = lambda wildcards, attempt: attempt * 3000,
        time = lambda wildcards, attempt: as_minutes(hours=2**(attempt - 1))
    shell:
        'scripts/find_signature_distances_from_ref.py {input} --output {output} --signature-type SourmashSignature'


rule stream_draff:
    input:
        'results/{sample_accession}/{producer}/to-draff-stream.pipe.fastq'
    output:
        stream = 'results/{sample_accession}/draff-stream/from-{producer}.sigs.gz',
        final = 'results/{sample_accession}/draff-stream/from-{producer}.final.sig'
    resources:
        mem = 8000,
        time = lambda _: as_minutes(hours=36)
    log:
        'logs/{sample_accession}/draff-stream/from-{producer}.log'
    shell:
        'goetia sketch draff -W {config[K]} -K 9 '
        '-S BitStorage -x 4e5 --interval {config[sketch_distance_interval]} '
        '--save-sig {output.final} --save-stream /dev/stdout --save-stream-tick 5 '
        '--pairing-mode single -i {input} '
        '2> {log} | gzip --rsyncable -c > {output.stream}'


rule stream_draff_distances:
    input:
        'results/{sample_accession}/draff-stream/from-{producer}.sigs.gz'
    output:
        'results/{sample_accession}/draff-stream/from-{producer}.stream-distances.csv'
    resources:
        mem = lambda wildcards, attempt: attempt * 3000,
        time = lambda wildcards, attempt: as_minutes(hours=2**(attempt - 1))
    shell:
        'scripts/find_signature_distances.py {input} --output {output} '
        '--signature-type DraffSignature --distance-metric cosine canberra '
        'correlation braycurtis'


rule stream_draff_ref_distances:
    input:
        'results/{sample_accession}/draff-stream/from-{producer}.sigs.gz'
        #ref = 'results/{sample_accession}/sourmash-stream/from-{producer}.final.sig'
    output:
        'results/{sample_accession}/draff-stream/from-{producer}.ref-distances.csv'
    resources:
        mem = lambda wildcards, attempt: attempt * 3000,
        time = lambda wildcards, attempt: as_minutes(hours=2**(attempt - 1))
    shell:
        'scripts/find_signature_distances_from_ref.py {input} --output {output} '
        '--signature-type DraffSignature --distance-metric cosine canberra '
        'correlation braycurtis'


rule stream_diginorm_draff:
    input:
        expand('results/{sample_accession}/{consumer}/from-diginorm-{diginorm_coverage}.{dtype}-distances.csv',
               sample_accession = samples.index.unique(),
               consumer = ['draff-stream'],
               diginorm_coverage = config['diginorm_coverage'],
               dtype = ['stream', 'ref'])


rule stream_diginorm_sourmash:
    input:
        expand('results/{sample_accession}/{consumer}/from-diginorm-{diginorm_coverage}.{dtype}-distances.csv',
               sample_accession = samples.index.unique(),
               consumer = ['sourmash-stream'],
               diginorm_coverage = config['diginorm_coverage'],
               dtype = ['stream', 'ref'])


rule stream_diginorm_signatures:
    input:
        expand('results/{sample_accession}/{consumer}/from-diginorm-{diginorm_coverage}.{dtype}-distances.csv',
               sample_accession = samples.index.unique(),
               consumer = ['sourmash-stream', 'draff-stream'],
               diginorm_coverage = config['diginorm_coverage'],
               dtype = ['stream', 'ref'])


rule stream_unfiltered_draff:
    input:
        expand('results/{sample_accession}/{consumer}/from-unfiltered.{dtype}-distances.csv',
               sample_accession = samples.index.unique(),
               consumer = ['draff-stream'],
               dtype = ['stream', 'ref'])


rule stream_unfiltered_sourmash:
    input:
        expand('results/{sample_accession}/{consumer}/from-unfiltered.{dtype}-distances.csv',
               sample_accession = samples.index.unique(),
               consumer = ['sourmash-stream'],
               dtype = ['stream', 'ref'])


rule stream_unfiltered_signatures:
    input:
        expand('results/{sample_accession}/{consumer}/from-unfiltered.{dtype}-distances.csv',
               sample_accession = samples.index.unique(),
               consumer = ['sourmash-stream', 'draff-stream'],
               dtype = ['stream', 'ref'])
