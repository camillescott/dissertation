#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : common.snakefile
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 08.04.2020


from datetime import timedelta
import pandas as pd


def read_sample_csv(path):
    return pd.read_csv(path, index_col='sample_accession', comment='#')


def build_sample_urls(samples):
    return {sample_acc: list(samples.loc[sample_acc].fastq_ftp) for sample_acc in samples.index.unique()}


def paired_stream_url_subst(accession, samples):
    url = samples.loc[accession]['fastq_ftp']
    if ';' in url:
        l, _, r = url.partition(';')
        return f'<(curl -sL {l.strip(";")}) <(curl -sL {r.strip(";")})'
    raise ValueError('Should be paired.')


def paired_url_split(accession, samples):
    url = samples.loc[accession]['fastq_ftp']
    if ';' in url:
        l, _, r = url.partition(';')
        return [l.strip(';'), r.strip(';')]
    raise ValueError('Should be paired.')


def get_sample_files(accession):
    return [os.path.join('data', 'runs', filename) for filename in samples.loc[accession].fastq_filename]


def as_minutes(**kwargs):
    return timedelta(**kwargs).seconds // 60

