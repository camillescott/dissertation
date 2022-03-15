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

from ficus import FigureManager, SubFigureManager
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

######################
#
# Data Generation
#
######################


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


rule stream_solid:
    '''
    Filter out reads that don't have the specified proportion of k-mers
    with count over the given threshold.
    '''
    conda: 'envs/goetia.yml'
    input:
        left  = 'data/fastx/{accession}.1.fq.gz',
        right = 'data/fastx/{accession}.2.fq.gz'
    output:
        pipe('results/chap1/solid-filter/{accession}.to-{consumer}.pipe.fastq')
    resources:
        mem       = 8000,
        time      = lambda _: as_minutes(hours=8),
    params:
        K         = 31,
        min_count = 3,
        min_prop  = 0.8
    threads: 1
    log: 'logs/solid-filter/{accession}.to-{consumer}.pipe.log'
    shell:
        'goetia filter solid -i {input.left} {input.right} --pairing-mode split '
        '-K {params.K} -C {params.min_count} -P {params.min_prop} '
        '-S ByteStorage -H CanLemireShifter -x 2e9 -N 4 '
        '--interval 1000000 -o {output} > {log} 2>&1'


rule solid_cdbg_build:
    conda: 'envs/goetia.yml'
    input:
        'results/chap1/solid-filter/{accession}.to-solid-cdbg.pipe.fastq'
    output:
        'results/chap1/solid-cdbg-build/{accession}/goetia.cdbg.components.json'
    threads: 3
    resources:
        mem = 32000,
        time = lambda _: as_minutes(hours=8)
    log: 'logs/solid-cdbg-build/{accession}/goetia.log'
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
        --results-dir results/chap1/solid-cdbg-build/{wildcards.accession}/ \
        --pairing-mode single -i {input}
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


rule all_solid_cdbg_build_txomes:
    input:
        expand('results/chap1/solid-cdbg-build/{accession}/goetia.cdbg.components.json',
               accession = TXOMIC_SAMPLES.index.unique())


rule all_hash_stream_baseline:
    input:
        expand('results/chap1/hash-stream-baseline/{accession}/metrics.json',
                accession = TXOMIC_SAMPLES.index.unique())


rule all_dbg_stream_baseline:
    input:
        expand('results/chap1/dbg-stream-baseline/{accession}/metrics.json',
                accession = TXOMIC_SAMPLES.index.unique())


######################
#
# Figure Generation
#
######################


@mpl.ticker.FuncFormatter
def numerize_fmtr(x, pos):
    return numerize.numerize(x)


def get_diffs(df):
    return pd.DataFrame({'d_t': df['elapsed_s'][1:].values - df['elapsed_s'][:-1].values,
                         'd_bytes': df['bytes_read'][1:].values - df['bytes_read'][:-1].values,
                         't': df['elapsed_s'][1:]})


def despine(*args, **kwargs):
    for ax in args:
        sns.despine(ax=ax, offset=10, trim=True, **kwargs)
        ax.yaxis.grid(ls='--')


def sample_desc(sample_name, metadata, n_reads=None, primary='species', multiline=False):
    sample_meta = metadata.loc[sample_name]
    species = ' '.join((f'$\it{{{w}}}$' for w in sample_meta.scientific_name.split()))
    n_reads = numerize.numerize(sample_meta.read_count) if n_reads is None else numerize.numerize(int(n_reads))
    div = '\n' if multiline else ''
    
    if primary == 'species':
        return f'{species} {div}({sample_name}, {n_reads} {sample_meta.library_strategy} {sample_meta.library_selection} reads)'
    else:
        return f'{sample_name} {div}({n_reads} {sample_meta.library_strategy} {sample_meta.library_selection} reads)'


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


def get_cdbg_comps_df():
    files = sorted(glob.glob('results/chap1/cdbg-build/*/goetia.cdbg.components.json'))
    
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
    return comps_df


def get_cdbg_build_metrics():

    files = sorted(glob.glob('results/chap1/cdbg-build/*/goetia.cdbg.stats.json'))
    samples = set([f.split('/')[4] for f in files])
    
    metrics_df = []
    for f in files:
        print(f)
        try:
            df = pd.read_json(f)
            df, prop_cols = normalize_metrics(df)
        except ValueError:
            pass
        else:
            #df['t_norm'] = df['t'] / df['t'].max()
            metrics_df.append(df)
    metrics_df = pd.concat(metrics_df).reset_index(drop=True)
    
    return metrics_df


def get_cdbg_frag_df():
    files = sorted(glob.glob('results/chap1/cdbg-build/*/goetia.cdbg.unitigs.bp.json'))
    
    frag_df = []
    for f in files:
        print(f)
        try:
            fdf = pd.read_json(f)
            fdf['t_norm'] = fdf['t'] / fdf['t'].max()
            frag_df.append(fdf)
        except ValueError as e:
            print(f'Error reading {f}: {e}')
            pass
    frag_df = pd.concat(frag_df).reset_index(drop=True)
    return frag_df


def convert_frag_to_proportion(input_df):
    bin_names = ['[31,50)', '[50,100)', '[100,200)', '[200,400)', '[400,800)', '[800,Inf)']
    output_df = input_df.copy()
    output_df[bin_names] = output_df[bin_names].div(output_df[bin_names].sum(axis=1), axis='rows')
    return output_df


def get_solid_cdbg_build_metrics():

    files = sorted(glob.glob('results/chap1/solid-cdbg-build/*/goetia.cdbg.stats.json'))
    samples = set([f.split('/')[4] for f in files])
    
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
    metrics_df['sample_name'] = metrics_df['sample_name'].str.partition('.')[0]
    
    return metrics_df


def get_solid_cdbg_comps_df():
    files = sorted(glob.glob('results/chap1/solid-cdbg-build/*/goetia.cdbg.components.json'))
    
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
    comps_df['sample_name'] = comps_df['sample_name'].str.partition('.')[0]
    return comps_df


def get_solid_cdbg_frag_df():
    files = sorted(glob.glob('results/chap1/solid-cdbg-build/*/goetia.cdbg.unitigs.bp.json'))
    
    frag_df = []
    for f in files:
        print(f)
        try:
            fdf = pd.read_json(f)
            fdf['t_norm'] = fdf['t'] / fdf['t'].max()
            frag_df.append(fdf)
        except ValueError as e:
            print(f'Error reading {f}: {e}')
            pass
    frag_df = pd.concat(frag_df).reset_index(drop=True)
    frag_df['sample_name'] = frag_df['sample_name'].str.partition('.')[0]
    return frag_df


rule chap_one_results_figure_one:
    output:
        'index/figure/chap1/chap1-results-dl-build-speed.png'
    params:
        n_samples = 8
    run:
        import matplotlib.pyplot as plt
        plt.switch_backend('Agg')
        import random

        files = sorted(glob.glob('results/chap1/SAM*.curl.txt'))
        samples = list(set([os.path.basename(f).split('.')[0] for f in files]))

        stream_df = []
        for accession in samples[:params.n_samples]:
            print(f'Reading {accession}...')
            left_fn = f'results/chap1/{accession}.1.curl.txt'
            right_fn = f'results/chap1/{accession}.2.curl.txt'
            ldf = get_diffs(pd.read_table(left_fn, delim_whitespace=True, names=['elapsed_s', 'bytes_read']))
            rdf = get_diffs(pd.read_table(right_fn, delim_whitespace=True, names=['elapsed_s', 'bytes_read']))
            
            bins = list(range(0, int(max(ldf.t.max(), rdf.t.max())), 20))
            
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
            #random.seed(42)
            #subsampled = random.sample(list(stream_df.accession.unique()), 8)
        
            sns.lineplot(data=stream_df, x='t', y='MiB/s', hue='accession', lw=1.5, legend=True, ax=ax)
            ax.set_title('Data rate of compaction stream piped from $\mathtt{curl}$')
            ax.set_ylabel('Stream Data Rate (MiB/s)')
            ax.set_xlabel('Time in Stream (seconds)')
            ax.legend(bbox_to_anchor=(1.05, .5), loc='center left', title='Sample Accession')
            despine(ax)
            #sns.despine(fig=fig, offset=10, trim=True)
            #ax.yaxis.grid(ls='--')


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
             FigureManager(filename=output[0], figsize=(12,8)) as (fig, ax):
            
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
                           nrows=len(files)//2,
                           ncols=2) as (fig, axs):
            
            for i, sample_name in enumerate(comps_df.sample_name.unique()):
                ax = axs[i // 2, i % 2]
            
                mdf = comps_df.query(f'sample_name == "{sample_name}"')
                g = sns.lineplot(data=mdf, x='t_norm', y='n_components', lw=2, color=sns.color_palette()[0], ax=ax)
                sax = g.axes.twinx()
                sns.lineplot(data=mdf, x='t_norm', y='max', lw=2, color=sns.color_palette()[1], ax=sax)
                ax.set_ylabel('Components')
                ax.set_xlabel('Normalized Position in Stream')
                ax.set_title(sample_name)
                sax.set_ylabel('Max Component Size')
                
                ax.yaxis.set_major_formatter(numerize_fmtr)
                sax.yaxis.set_major_formatter(numerize_fmtr)

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
             FigureManager(filename=output[0], ncols=2, figsize=(12,6)) as (fig, axs):
            
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
            axs[0].set_ylabel('$k$-mers/s')
            
            sns.despine(ax=axs[1], offset=10, trim=True)
            axs[1].yaxis.grid(ls='--')
            axs[1].yaxis.set_major_formatter(numerize_fmtr)
            axs[1].set_ylabel('sequences/s')


rule chap_one_results_figure_five:
    output:
        'index/figure/chap1/chap1-results-solid-human.png'
    run:
        import matplotlib.pyplot as plt
        plt.switch_backend('Agg')
        from matplotlib.lines import Line2D

        raw_cdbg_metrics_df = get_cdbg_build_metrics()
        raw_comps_df = get_cdbg_comps_df()

        solid_cdbg_metrics_df = get_solid_cdbg_build_metrics()
        solid_comps_df = get_solid_cdbg_comps_df()

        raw_cdbg_metrics_df['prefilter'] = 'raw'
        raw_comps_df['prefilter'] = 'raw'

        solid_cdbg_metrics_df['prefilter'] = 'solid'
        solid_comps_df['prefilter'] = 'solid'

        prefilter_metrics_df = pd.concat((raw_cdbg_metrics_df, solid_cdbg_metrics_df)).reset_index(drop=True)
        prefilter_comps_df = pd.concat((raw_comps_df, solid_comps_df)).reset_index(drop=True)

        samples = TXOMIC_SAMPLES.query('scientific_name == "Homo sapiens"').index

        with sns.axes_style("ticks"), \
             SubFigureManager(filename=output[0],
                              figsize=(12,len(samples) * 4),  
                              subfigs='row',
                              nrows=len(samples), 
                              ncols=3) as (fig, subfigs, subaxes):
            
            for row, (sample_name, subfig, subax) in enumerate(zip(samples, subfigs, subaxes)):
                lax, cax, rax = subax
                df = prefilter_comps_df.query(f'sample_name == "{sample_name}"')
                mdf = prefilter_metrics_df.query(f'sample_name == "{sample_name}"')

                subfig.suptitle(sample_desc(sample_name, TXOMIC_SAMPLES, n_reads=mdf['seq_t'].max()))
                
                sns.lineplot(data=df, x='t_norm', y='n_components', hue='prefilter', legend=False, lw=2, ax=lax)
                lax.set_ylabel('Components')
                lax.set_xlabel('Normalized Position in Stream')

                sns.lineplot(data=df, x='t_norm', y='max', hue='prefilter', legend=False, lw=2, ax=cax)
                cax.set_ylabel('Max Component Size')
                cax.set_xlabel('Normalized Position in Stream')
                
                sns.lineplot(data=mdf, x='t_norm', y='n_unique_kmers', hue='prefilter', lw=2, ax=rax)
                rax.set_ylabel('Unique $k$-mers')
                rax.set_xlabel('Normalized Position in Stream')
                
                rax.legend(bbox_to_anchor=(1.05, .5), loc='center left', title='Sequence Prefilter')

                despine(lax, cax, rax)
                lax.yaxis.set_major_formatter(numerize_fmtr)
                cax.yaxis.set_major_formatter(numerize_fmtr)
                rax.yaxis.set_major_formatter(numerize_fmtr)


rule chap_one_results_figure_six:
    output:
        'index/figure/chap1/chap1-results-human-frag-and-metrics.png'
    run:
        import matplotlib.pyplot as plt
        plt.switch_backend('Agg')
        from matplotlib.lines import Line2D

        metrics_df = get_cdbg_build_metrics()
        frag_df = get_cdbg_frag_df()
        prop_frag_df = convert_frag_to_proportion(frag_df).melt(id_vars=['sample_name', 't', 't_norm'], var_name='Unitig Length Bin')
        frag_df = frag_df.melt(id_vars=['sample_name', 't', 't_norm'], var_name='Unitig Length Bin')

        human  = TXOMIC_SAMPLES.query('scientific_name == "Homo sapiens"').index

        with sns.axes_style("ticks"), \
             SubFigureManager(filename=output[0],
                              subfigs='col',
                              ncols=len(human),
                              nrows=2,
                              subplot_kwds=dict(sharex=True),
                              figsize=(12,8)) as (fig, subfigs, subaxes):
            
            for i, (sample_name, subfig, subax) in enumerate(zip(human, subfigs, subaxes)):
                tax, bax = subax
                legend = True if i == len(human) - 1 else False
                
                fdf = prop_frag_df.query(f'sample_name == "{sample_name}"')
                mdf = metrics_df.query(f'sample_name == "{sample_name}"')[['t', 't_norm', 'seq_t', 'n_full', 'n_tips', 'n_islands', 'n_dnodes']]
                n_reads = mdf['seq_t'].max()
                mdf = mdf.melt(id_vars=['t', 't_norm', 'seq_t'], var_name='Node Type', value_name='Unitig Count')
                
                with sns.color_palette("viridis"):
                    sns.lineplot(data=fdf, x='t_norm', y='value',  hue='Unitig Length Bin', legend=legend, lw=1.5, ax=tax)
                
                sns.lineplot(data=mdf, x='t_norm', y='Unitig Count',  hue='Node Type', legend=legend, lw=1.5, ax=bax)
                bax.set_xlabel('Normalized Position in Stream')
                
                if i == 0:
                    tax.set_ylabel('Proportion of Unitigs in Bin')
                else:
                    tax.set_ylabel('')
                    bax.set_ylabel('')
                    tax.get_yaxis().set_visible(False)
                    
                tax.set_ylim(ymax=1.0)
                bax.yaxis.set_major_formatter(numerize_fmtr)
                
                subfig.suptitle(sample_desc(sample_name, TXOMIC_SAMPLES, primary='accession', n_reads=n_reads, multiline=True))
            
                if legend:
                    tax.legend(loc='upper right', title='Unitig Length Bin')
                    bax.legend(['Full', 'Tip', 'Island', 'Decision'])
               

rule chap_one_results_figure_seven:
    output:
        'index/figure/chap1/chap1-results-solid-human-frag.png'
    run:
        import matplotlib.pyplot as plt
        plt.switch_backend('Agg')
        from matplotlib.lines import Line2D

        def melt_frag_df(df):
            return df.melt(id_vars=['sample_name', 't', 't_norm'], var_name='Unitig Length Bin')

        raw_frag_df = convert_frag_to_proportion(get_cdbg_frag_df())
        solid_frag_df = convert_frag_to_proportion(get_solid_cdbg_frag_df())
        metrics_df = get_cdbg_build_metrics()

        samples = TXOMIC_SAMPLES.query('scientific_name == "Homo sapiens"').index
        bin_names = ['[31,50)', '[50,100)', '[100,200)', '[200,400)', '[400,800)', '[800,Inf)']
        cuts = np.linspace(0, 1, 20)

        with sns.axes_style("ticks"), \
             SubFigureManager(filename=output[0],
                              figsize=(12,len(samples) * 4),  
                              subfigs='row',
                              nrows=len(samples),
                              ncols=3) as (fig, subfigs, subaxes):
            
            from matplotlib.lines import Line2D
            
            for row, (sample_name, subfig, subax) in enumerate(zip(samples, subfigs, subaxes)):
                lax, cax, rax = subax
                mdf = metrics_df.query(f'sample_name == "{sample_name}"')
                
                solid_df = solid_frag_df.query(f'sample_name == "{sample_name}"')
                norm_solid_df = solid_df.groupby(pd.cut(solid_df.t_norm, labels=cuts[1:], bins=cuts)).mean()
                
                raw_df = raw_frag_df.query(f'sample_name == "{sample_name}"')
                norm_raw_df = raw_df.groupby(pd.cut(raw_df.t_norm, labels=cuts[1:], bins=cuts)).mean()
                
                diff_df = (norm_solid_df[bin_names] - norm_raw_df[bin_names]).reset_index()
                
                subfig.suptitle(sample_desc(sample_name, TXOMIC_SAMPLES, n_reads=mdf['seq_t'].max()))
                
                with sns.color_palette("hls", 7):
                    sns.lineplot(data=melt_frag_df(raw_df),
                                 x='t_norm',
                                 y='value', 
                                 hue='Unitig Length Bin', legend=False, lw=1.5, ax=lax)

                    sns.lineplot(data=melt_frag_df(solid_df),
                                 x='t_norm',
                                 y='value', 
                                 hue='Unitig Length Bin', legend=False, lw=1.5, ax=cax)

                    sns.lineplot(data=diff_df.melt(id_vars=['t_norm'], var_name='Unitig Length Bin'),
                                 x='t_norm',
                                 y='value', 
                                 hue='Unitig Length Bin', legend=True, lw=1.5, ax=rax)

                lax.set_ylim(ymin=0, ymax=1.0)
                cax.set_ylim(ymin=0, ymax=1.0)
                rax.set_ylim(ymin=-.3, ymax=.3)
                
                rax.legend(bbox_to_anchor=(1.05, .5), loc='center left', title='Unitig Length Bin')

                lax.yaxis.grid(ls='--')
                cax.yaxis.grid(ls='--')
                rax.yaxis.grid(ls='--')
                
                lax.set_title('Raw')
                cax.set_title('Solid')
                rax.set_title('Difference (Solid - Raw)')
                
                lax.set_xlabel('Normalized Position in Stream')
                cax.set_xlabel('Normalized Position in Stream')
                rax.set_xlabel('Normalized Position in Stream')


rule all_figures:
    input:
        rules.chap_one_results_figure_one.output,
        rules.chap_one_results_figure_two.output,
        rules.chap_one_results_figure_three.output,
        rules.chap_one_results_figure_four.output,
        rules.chap_one_results_figure_five.output,
        rules.chap_one_results_figure_six.output,
        rules.chap_one_results_figure_seven.output
