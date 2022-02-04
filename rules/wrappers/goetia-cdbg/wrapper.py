__author__ = "Camille Scott"
__copyright__ = "Copyright 2022, Camille Scott"
__email__ = "cswel@ucdavis.edu"
__license__ = "MIT"


import os
from snakemake.shell import shell

opts = []
opts.append('--results-dir {0}'.format(
            snakemake.params.get('results_dir',
                                 snakemake.rule + '.cdbg.results')))

if snakemake.params.get('track_cdbg_metrics'):
    opts.append('--track-cdbg-metrics')
if snakemake.params.get('track_cdbg_components'):
    opts.append('--track-cdbg-components')
    opts.append('--component-sample-size {0}'.format(
                snakemake.params.get('sample_size', 10000)))

if snakemake.params.get('track_unitig_bp'):
    opts.append('--track-unitig-bp')

normalize = snakemake.params.get('normalize', False)
if normalize:
    opts.append('--normalize')
    if isinstance(normalize, int) and not isinstance(normalize, bool):
        opts.append(str(normalize))

save_cdbg = snakemake.params.get('save_cdbg', False)
if save_cdbg:
    opts.append('--save-cdbg')
    if isinstance(save_cdbg, str):
        opts.append(save_cdbg)
    if snakemake.params.get('save_cdbg_format'):
        opts.append('--save-cdbg-format')
        opts.extend(snakemake.params.get('save_cdbg_format'))

if snakemake.params.get('validate'):
    opts.append('--validate')

opts.append('--storage {0}'.format(
            snakemake.params.get('storage_type', 'PHMapStorage')))
opts.append('--interval {0}'.format(
            snakemake.params.get('interval', 5000000)))

opts.append('-K {0}'.format(
            snakemake.params.get('ksize', 31)))

if snakemake.params.get('extra'):
    opts.append(snakemake.params.get('extra'))


if not hasattr(snakemake.input, 'r1'):
    # assume single
    opts.append('--pairing-mode single')
    if isinstance(snakemake.input, str):
        inputs = snakemake.input
    else:
        inputs = ' '.join(snakemake.input)
else:
    if not hasattr(snakemake.input, 'r2'):
        raise ValueError('r1 requires r2 to also be defined as input.')
    opts.append('--pairing-mode split')
    if isinstance(snakemake.input.r1, str):
        inputs = '{0} {1}'.format(snakemake.input.r1,
                                  snakemake.input.r2)
    else:
        r1 = list(snakemake.input.r1)
        r2 = list(snakemake.input.r2)
        files = []
        while r2:
            files.append(r1.pop())
            files.append(r2.pop())
        inputs = ' '.join(files)

print(opts)
opts = ' '.join(opts)
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell('goetia cdbg build -H FwdLemireShifter {opts} '
      '-i {inputs} {log}')
