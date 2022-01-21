#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : common.snakefile
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 08.04.2020


from datetime import timedelta
import pandas as pd


samples = pd.read_csv('samples.csv', index_col='sample_accession')
urls    = {sample_acc: list(samples.loc[sample_acc].fastq_ftp) for sample_acc in samples.index.unique()}

# Subset the accessions to use
small_samples = (samples.groupby('sample_accession')['read_count'].sum() / 2 < 250000000)
samples = samples[small_samples]
samples = pd.concat([samples.query('library_source == "TRANSCRIPTOMIC"').head(config['n_txomes'] * 2), 
                     samples.query('library_source == "GENOMIC"').head(config['n_genomes'] * 2)])


def paired_stream_url_subst(accession):
    url = samples.loc[accession]['fastq_ftp']
    if ';' in url:
        l, _, r = url.partition(';')
        return f'<(curl -sL {l.strip(";")}) <(curl -sL {r.strip(";")})'
    raise ValueError('Should be paired.')


def paired_url_split(accession):
    url = samples.loc[accession]['fastq_ftp']
    if ';' in url:
        l, _, r = url.partition(';')
        return [l.strip(';'), r.strip(';')]
    raise ValueError('Should be paired.')


def get_sample_files(accession):
    return [os.path.join('data', 'runs', filename) for filename in samples.loc[accession].fastq_filename]


def as_minutes(**kwargs):
    return timedelta(**kwargs).seconds // 60

