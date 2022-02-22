#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2022
# File   : common.snakefile
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 01.02.2022

import glob
import os

import matplotlib.pyplot as plt
plt.switch_backend('Agg')

from ficus import FigureManager
import ijson
import matplotlib as mpl
import numpy as np
import pandas as pd
import seaborn as sns

from numerize import numerize


GENOMIC_SAMPLES = read_sample_csv('config/chap1-genomes.csv')
GENOMIC_URLS = build_sample_urls(GENOMIC_SAMPLES)

TXOMIC_SAMPLES = read_sample_csv('config/chap1-txomes.csv')
TXOMIC_URLS = build_sample_urls(TXOMIC_SAMPLES)

#print('Chapter One Txomic Accessions: ', ' '.join(TXOMIC_SAMPLES.index.unique()))


@mpl.ticker.FuncFormatter
def numerize_fmtr(x, pos):
    return numerize.numerize(x)


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


def parse_components_metrics(fp):
    backend = ijson.get_backend('yajl2')
    metrics = []
    #sizes = []
    for row in backend.items(fp, 'item'):
        metrics.append({k: row[k] for k in set(row.keys()) - set(('sizes',))})
        #sizes.append(sorted(row['sizes']))
    metrics_df = pd.DataFrame(metrics)
    #sizes_df = pd.DataFrame(sizes, index=metrics_df.index)
    return metrics_df


rule download_stream_cdbg_build:
    conda: 'envs/goetia.yml'
    input:
        left  = 'data/streams/{accession}.pipe.1.fq.gz',
        right = 'data/streams/{accession}.pipe.2.fq.gz'
    output:
        'results/chap1/cdbg-stream/{accession}/goetia.cdbg.stats.json'
    threads: 3
    resources:
        mem = 96000,
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


rule predownloaded_cdbg_build:
    conda: 'envs/goetia.yml'
    input:
        left  = 'data/fastx/{accession}.1.fq.gz',
        right = 'data/fastx/{accession}.2.fq.gz'
    output:
        'results/chap1/cdbg-build/{accession}/goetia.cdbg.components.json'
    threads: 3
    resources:
        mem = 32000,
        time = lambda _: as_minutes(hours=8)
    log: 'logs/cdbg-build/build-{accession}.log'
    params:
        storage_type     = 'PHMapStorage',
        interval         = 1000000,
        K                = 31,
        sample_size      = 10000,
        bins             = ' '.join((str(b) for b in [31, 50, 100, 200, 400, 800]))
    shell: '''
        goetia cdbg build -K {params.K} -S PHMapStorage -H FwdLemireShifter --interval {params.interval} \
        --track-cdbg-metrics --verbose \
        --track-cdbg-components --component-sample-size {params.sample_size} --cdbg-components-tick exp \
        --track-unitig-bp --unitig-bp-bins {params.bins} --unitig-bp-tick 10 \
        --results-dir results/chap1/cdbg-build/{wildcards.accession}/ \
        --pairing-mode split -i {input.left} {input.right}
    '''


rule hash_stream_baseline:
    conda: 'envs/goetia.yml'
    input:
        left  = 'data/fastx/{accession}.1.fq.gz',
        right = 'data/fastx/{accession}.2.fq.gz'
    output:
        'results/chap1/hash-stream-baseline/{accession}/metrics.json'
    threads: 3
    resources:
        mem = 4000,
        time = lambda _: as_minutes(hours=1)
    params:
        K = 31,
        interval = 1000000
    shell: '''
        goetia utils hash-stream -K {params.K} --interval {params.interval} \
        --metrics {output} \
        --pairing-mode split -i {input.left} {input.right}
    '''

rule dbg_stream_baseline:
    conda: 'envs/goetia.yml'
    input:
        left  = 'data/fastx/{accession}.1.fq.gz',
        right = 'data/fastx/{accession}.2.fq.gz'
    output:
        'results/chap1/dbg-stream-baseline/{accession}/metrics.json'
    threads: 3
    resources:
        mem = 4000,
        time = lambda _: as_minutes(hours=1)
    params:
        K = 31,
        interval = 1000000
    shell: '''
        goetia utils hash-stream -K {params.K} --dbg --interval {params.interval} \
        --metrics {output} \
        --pairing-mode split -i {input.left} {input.right}
    '''


rule all_download_stream_cdbg_build_txomes:
    input:
        expand('results/chap1/cdbg-stream/{accession}/goetia.cdbg.stats.json',
               accession = TXOMIC_SAMPLES.index.unique())


rule all_predownloaded_cdbg_build_txomes:
    input:
        expand('results/chap1/cdbg-build/{accession}/goetia.cdbg.components.json',
               accession = TXOMIC_SAMPLES.index.unique())


rule all_hash_stream_baseline:
    input:
        expand('results/chap1/hash-stream-baseline/{accession}/metrics.json',
                accession = TXOMIC_SAMPLES.index.unique())


rule all_dbg_stream_baseline:
    input:
        expand('results/chap1/dbg-stream-baseline/{accession}/metrics.json',
                accession = TXOMIC_SAMPLES.index.unique())


rule chap_one_results_figure_one:
    output:
        'index/figure/chap1/chap1-results-dl-build-speed.png'
    run:
        import matplotlib.pyplot as plt
        plt.switch_backend('Agg')
        import random

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
            random.seed(42)
            subsampled = random.sample(list(stream_df.accession.unique()), 8)
        
            sns.lineplot(data=stream_df, x='t', y='MiB/s', hue='accession', lw=1, legend=True, ax=ax)
            ax.set_title('Data rate of compaction stream piped from $\mathtt{curl}$')
            ax.set_ylabel('Stream Data Rate (MiB/s)')
            ax.set_xlabel('Time in Stream (seconds)')
            sns.despine(fig=fig, offset=10, trim=True)
            ax.yaxis.grid(ls='--')


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
            sns.despine(fig=fig, offset=10, trim=True)
            ax.yaxis.grid(ls='--')


rule chap_one_results_figure_three:
    output:
        'index/figure/chap1/chap1-results-txome-comps.png'
    params:
        n_subplots = 8
    run:
        import matplotlib.pyplot as plt
        plt.switch_backend('Agg')
        from matplotlib.lines import Line2D

        files = sorted(glob.glob('results/chap1/cdbg-build/*/goetia.cdbg.components.json'))
        files = files[:params.n_subplots]
        comps_df = []
        for f in files:
            print(f)
            try:
                with open(f) as fp:
                    mdf = parse_components_metrics(fp)
                    mdf['t_norm'] = mdf['t'] / mdf['t'].max()
                    comps_df.append(mdf)
            except ijson.IncompleteJSONError:
                print(f'JSON error in {f}')
                pass
        comps_df = pd.concat(comps_df).reset_index(drop=True)
        comps_df['log_n_components'] = np.log(comps_df['n_components'])

        with sns.axes_style("ticks"), \
             FigureManager(filename=output[0],
                           figsize=(8,12),
                           tight_layout=True, 
                           nrows=len(files)//2,
                           ncols=2) as (fig, axs):
            
            for i, sample_name in enumerate(comps_df.sample_name.unique()):
                ax = axs[i // 2, i % 2]
            
                #mdf = comps_df.query('sample_name == "SAMN09758735"')
                mdf = comps_df.query(f'sample_name == "{sample_name}"')
                #df = mdf.reset_index().melt(value_vars=['min', 'max', 'n_components'], id_vars=['sample_name', 't'])
                g = sns.lineplot(data=mdf, x='t_norm', y='n_components', lw=2, color=sns.color_palette()[0], ax=ax)
                sax = g.axes.twinx()
                sns.lineplot(data=mdf, x='t_norm', y='max', lw=2, color=sns.color_palette()[1], ax=sax)
                ax.set_ylabel('Components')
                ax.set_xlabel('Normalized Position in Stream')
                ax.set_title(sample_name)
                sax.set_ylabel('Max Component Size')

                g.legend(handles=[Line2D([], [], marker='_', color=sns.color_palette()[0], label='Components'), 
                                  Line2D([], [], marker='_', color=sns.color_palette()[1], label='Max Component Size')],
                        loc='lower right')
                sns.despine(ax=ax, offset=10, trim=True, right=False)
                sns.despine(ax=sax, offset=10, trim=True, right=False)
                ax.yaxis.grid(ls='--')


rule chap_one_results_figure_four:
    output:
        'index/figure/chap1/chap1-results-rate-comp.png'
    run:
        def calc_rates(df):
            return pd.DataFrame({'kmers/s': (df['t'].values[1:] - df['t'].values[:-1]) / df['rt_elapsed_interval'].values[1:],
                                 'seqs/s': (df['seq_t'].values[1:] - df['seq_t'].values[:-1]) / df['rt_elapsed_interval'].values[1:],
                                 't': df['t'].values[1:],
                                 'sample_name': df['sample_name'].values[1:],
                                 't_norm': df['t_norm'].values[1:],
                                 'rt_elapsed_total': df['rt_elapsed_total'][1:]})

        def read_metrics_files(files):
            metrics_df = []
            for f in files:
                print(f)
                try:
                    df = pd.read_json(f)
                except ValueError as e:
                    print(f'Error reading {f}: {e}')
                    pass
                else:
                    df['t_norm'] = df['t'] / df['t'].max()
                    metrics_df.append(calc_rates(df))
            metrics_df = pd.concat(metrics_df).reset_index(drop=True)
            return metrics_df

        cdbg_metrics_files =  sorted(glob.glob('results/chap1/cdbg-build/*/goetia.cdbg.stats.json'))
        dbg_baseline_files = sorted(glob.glob('results/chap1/dbg-stream-baseline/*/metrics.json'))
        hash_baseline_files = sorted(glob.glob('results/chap1/hash-stream-baseline/*/metrics.json'))

        cdbg_metrics = read_metrics_files(cdbg_metrics_files)
        dbg_metrics = read_metrics_files(dbg_baseline_files)
        hash_metrics = read_metrics_files(hash_baseline_files)

        hash_metrics['operation'] = 'hash'
        dbg_metrics['operation'] = 'dbg'
        cdbg_metrics['operation'] = 'cdbg'

        rates_df = pd.concat((hash_metrics, dbg_metrics, cdbg_metrics))

        with sns.axes_style("ticks"), \
             FigureManager(filename=output[0], tight_layout=True, ncols=2, figsize=(12,6)) as (fig, axs):
            
            mean_rates = data=rates_df.groupby(['sample_name', 'operation']).mean().reset_index()
            sns.boxplot(x='operation', y='kmers/s', data=mean_rates, ax=axs[0])
            sns.swarmplot(x='operation', y='kmers/s', data=mean_rates, color=".25", ax=axs[0])
            
            sns.boxplot(x='operation', y='seqs/s', data=mean_rates, ax=axs[1])
            sns.swarmplot(x='operation', y='seqs/s', data=mean_rates, color=".25", ax=axs[1])
            
            #ax.legend(bbox_to_anchor=(1.05, .5), loc='center left', title='Sample Accession')
            #ax.set_title('Number of Components in cDBG')
            
            sns.despine(ax=axs[0], offset=10, trim=True)
            axs[0].yaxis.grid(ls='--')
            axs[0].yaxis.set_major_formatter(numerize_fmtr)
            
            sns.despine(ax=axs[1], offset=10, trim=True)
            axs[1].yaxis.grid(ls='--')
            axs[1].yaxis.set_major_formatter(numerize_fmtr)
