#!/usr/bin/env python3
"""
Submit this clustering script for sbatch to snakemake with:
    snakemake -j 99 --cluster slurm_scheduler.py
"""

import argparse
import os
import subprocess
import sys
import warnings

from snakemake.utils import read_job_properties


def make_dir(directory):
    """Make directory unless existing. Ignore error in the latter case."""
    try:
        os.makedirs(directory)
    except:
        pass


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("dependencies", nargs="*",
                        help="{{dependencies}} string given by snakemake\n")
    parser.add_argument("jobscript",
                        help="Snakemake generated shell script with commands to execute snakemake rule\n")
    args = parser.parse_args()

    if args.dependencies:
        dependencies = '-d afterok:' + ':'.join(args.dependencies)
    else:
        dependencies = ''

    job_properties = read_job_properties(args.jobscript)

    cluster_param = {}
    job_resources = job_properties["resources"]

    if not "mem" in job_resources:
        warnings.warn("Rule {rule} has no memory specified, set to default.".format(**job_properties))
   
    # do something useful with the threads
    cluster_param["threads"] = job_properties.get("threads",1)
    cluster_param['days']    = job_resources.get("days", job_properties["cluster"]['days'])
    cluster_param['hours']   = job_resources.get("hours", job_properties["cluster"]['hours'])
    cluster_param['mem']     = int(job_resources.get("mem", 1)) + 1 #GB + overhead

    cluster_param['name']    = job_properties['cluster']['name']
    cluster_param['account'] = job_properties['cluster']['account']
    cluster_param['partition'] = job_properties['cluster']['partition']
    cluster_param['email']   = job_properties['cluster']['email']

    logdir = job_properties['cluster']['logdir']
    make_dir(logdir)

    sbatch_cmd = "sbatch -A {account} -p {partition} --parsable -c {threads} --export=ALL "\
                 "{dependencies} "\
                 "-o {logdir}/slurm-%j.out "\
                 "--time={days:d}-{hours:02d}:00:00 --mem={mem}g "\
                 "--mail-type=FAIL,BEGIN,END --mail-user {email} "\
                 "--job-name={name} {script}".format(script=args.jobscript,
                                                     dependencies=dependencies,
                                                     logdir=logdir,
                                                     **cluster_param)

    print(sbatch_cmd, file=sys.stderr)
    popenrv = subprocess.Popen(sbatch_cmd,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT,
                               shell=True).communicate()

    # Snakemake expects only id of submitted job on stdout for scheduling
    # with {dependencies}
    try:
        print("%i" % int(popenrv[0].split()[-1]))
    except ValueError:
        print("Not a submitted job: %s" % popenrv[0])
        sys.exit(2)
