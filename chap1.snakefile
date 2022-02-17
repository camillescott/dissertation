#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2022
# File   : common.snakefile
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 01.02.2022

import glob
import os

from ficus import FigureManager
import numpy as np
import pandas as pd
import seaborn as sns


GENOMIC_SAMPLES = read_sample_csv('config/chap1-genomes.csv')
GENOMIC_URLS = build_sample_urls(GENOMIC_SAMPLES)

TXOMIC_SAMPLES = read_sample_csv('config/chap1-txomes.csv')
TXOMIC_URLS = build_sample_urls(TXOMIC_SAMPLES)


def get_diffs(df):
    return pd.DataFrame({'d_t': df['elapsed_s'][1:].values - df['elapsed_s'][:-1].values,
                         'd_bytes': df['bytes_read'][1:].values - df['bytes_read'][:-1].values,
                         't': df['elapsed_s'][1:]})


def normalize_metrics(data):
    data['t_norm'] = data['t'] / data['t'].max()

    cdbg_cols =['n_full', 'n_tips', 'n_islands', 'n_trivial', 'n_circular']
    for col in cdbg_cols:
        data[col + '_p'] = data[col] / data['n_unodes']
    data['dnode_p'] = data['n_dnodes'] / (data['n_unodes'] + data['n_dnodes'])
    data['kmer_p'] = data['n_unique_kmers'] / data['n_unique_kmers'].max()
    prop_cols = [col + '_p' for col in cdbg_cols] + [ 'kmer_p']
    
    return data, prop_cols


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


rule all_download_stream_cdbg_build:
    input:
        expand('results/chap1/cdbg-stream/{accession}/goetia.cdbg.stats.json',
               accession = TXOMIC_SAMPLES.index.unique())


rule chap_one_results_figure_one:
    output:
        'index/figure/chap1/chap1-results-dl-build-speed.png'
    run:
        import matplotlib.pyplot as plt
        plt.switch_backend('Agg')

        files = sorted(glob.glob('results/chap1/SAM*.curl.txt'))
        samples = set([os.path.basename(f).split('.')[0] for f in files])

        stream_df = []
        for accession in samples:
            print(f'Reading {accession}...')
            left_fn = f'results/chap1/{accession}.1.curl.txt'
            right_fn = f'results/chap1/{accession}.2.curl.txt'
            ldf = get_diffs(pd.read_table(left_fn, delim_whitespace=True, names=['elapsed_s', 'bytes_read']))
            rdf = get_diffs(pd.read_table(right_fn, delim_whitespace=True, names=['elapsed_s', 'bytes_read']))
            
            bins = list(range(0, int(max(ldf.t.max(), rdf.t.max())), 5))
            
            ldf = ldf.groupby(pd.cut(ldf.t, labels=bins[1:], bins=bins)).mean()
            ldf['bytes/s'] = ldf.d_bytes / ldf.d_t
            rdf = rdf.groupby(pd.cut(rdf.t, labels=bins[1:], bins=bins)).mean()
            rdf['bytes/s'] = rdf.d_bytes / rdf.d_t
            
            df = pd.DataFrame({'t': ldf.index,
                               'bytes/s': ldf['bytes/s'] + rdf['bytes/s'],
                               'accession': accession})
            df['MiB/s'] = df['bytes/s'] / 1048576

            stream_df.append(df.copy())
        stream_df = pd.concat(stream_df).reset_index(drop=True)

        with sns.axes_style("ticks"), FigureManager(filename=output[0], figsize=(12,8)) as (fig, ax):
        
            sns.lineplot(data=stream_df, x='t', y='MiB/s', hue='accession', lw=1, legend=True, ax=ax)
            ax.set_title('Data rate of compaction stream piped from $\mathtt{curl}$')
            ax.set_ylabel('Stream Data Rate (MiB/s)')
            ax.set_xlabel('Time in Stream (seconds)')


rule chap_one_results_figure_two:
    output:
        'index/figure/chap1/chap1-results-dnode-prop.png'
    run:
        import matplotlib.pyplot as plt
        plt.switch_backend('Agg')

        files = sorted(glob.glob('results/chap1/cdbg-stream/*/goetia.cdbg.stats.json'))

        metrics_df = []
        for f in files:
            print(f)
            try:
                df = pd.read_json(f)
                df, prop_cols = normalize_metrics(df)
            except ValueError:
                pass
            else:
                df['t_norm'] = df['t'] / df['t'].max()
                metrics_df.append(df)
        metrics_df = pd.concat(metrics_df).reset_index(drop=True)
        metrics_df['sample_name'] = metrics_df.sample_name.str.rpartition('.')[0]

        with sns.axes_style("ticks"), \
             FigureManager(filename=output[0], tight_layout=True, figsize=(12,8)) as (fig, ax):
            
            sns.lineplot(data=metrics_df, x='t_norm', y='dnode_p', hue='sample_name', lw=1, ax=ax)
            ax.set_ylabel('Decision Node Proportion')
            ax.set_xlabel('Normalized Position in Stream')
            ax.set_title('Proportion of cDBG nodes that are decision nodes')
            ax.legend(bbox_to_anchor=(1.05, .5), loc='center left', title='Sample Accession')
