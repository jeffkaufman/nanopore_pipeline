
#!/usr/bin/python3

import sys
import json
import glob
import os
import subprocess
import shutil
import argparse

from command_line_util import *

# equivalent to "conda activate basecall"
os.environ["PATH"] = os.pathsep.join([
  "/home/ec2-user/anaconda3/envs/basecall/bin", os.environ["PATH"]])

def select_basecall_model(report):
  basecall_models = json.loads(
    report['protocol_run_info']['meta_info']['tags']
    ['available basecall models']['array_value'])

  for model in basecall_models:
    if model.endswith('_sup.cfg'):
      return model

  die("Can't identify sup model from %r" % basecall_models)

def upload_public(fname):
  if '/' in fname:
    die('%s must be a local path' % fname)
  run(['aws', 's3', 'cp', fname, 's3://ontseq-res-public'])
  run(['aws', 's3api', 'put-object-acl', '--bucket', 'ontseq-res-public',
       '--key', fname,
       '--grant-read',
       'uri=http://acs.amazonaws.com/groups/global/AllUsers'])
  
def call_bases(run_dir, results_dir):
  report_fname, = glob.glob(os.path.join(run_dir, 'report_*.json'))
  with open(report_fname) as inf:
    report = json.load(inf)

  basecall_model = select_basecall_model(report)

  run([
    os.path.expanduser('~/ont-guppy/bin/guppy_basecaller'),
    '-i', os.path.join(run_dir, 'fast5/'),
    '-s', results_dir,
    '-c', os.path.join('ont-guppy/data', basecall_model),
    '--gpu_runners_per_device', '1',
    '--num_callers', '1',
    '--cpu_threads_per_caller', '8',
    '--device', 'cuda:0',
    '--verbose',
    '--barcode_kits', 'bw-4',
    '--disable_pings',
    '--allow_inferior_barcodes',
    '--require_barcodes_both_ends',
    '--num_barcoding_buffers', '128',
    '--arrangements_file', 'bw-3.cfg',
  ])

def start():
  parser = argparse.ArgumentParser(
    description='Call bases and determine consensus sequences')
  parser.add_argument('--upload-fastq', action='store_true')
  parser.add_argument('jobname')
  parser.add_argument('bararr')
  args = parser.parse_args()

  validate_jobname(args.jobname)

  tarfile = tarfile_name(args.jobname)
  if os.path.exists(tarfile):
    info('tarfile exists; skipping download from s3')
  else:
    run(['aws', 's3', 'cp',  s3_url(args.jobname), '.'])

  if not os.path.exists(args.bararr):
    die('expected %r to exist' % args.bararr)

  raw_dir = args.jobname
  if not os.path.exists(raw_dir):
    os.mkdir(raw_dir)
    run(['tar', '-xvf', tarfile, '-C', raw_dir])

  results_dir = "%sResults" % args.jobname
  if os.path.exists(results_dir):
    info('basecalling results already exist, skipping')
  else:
    os.mkdir(results_dir)

    report_fname, = glob.glob(os.path.join(raw_dir, '**/report_*.json'),
                              recursive=True)
    run_dir = os.path.dirname(report_fname)

    call_bases(run_dir, results_dir)

  bararr_in_results = os.path.join(results_dir, 'bararr.csv')
  if os.path.exists(bararr_in_results):
    info('bararr.csv already in output, skipping')
  else:
    shutil.copyfile(args.bararr, bararr_in_results)
    run([
      os.path.join(THISDIR, 'bc_csv_to_fasta.py'),
      results_dir,
    ])  # creates <results_dir>/samples.csv

  os.chdir(results_dir)

  if args.upload_fastq:
    if not os.path.exists('pass') and not os.path.exists('fail'):
      die('pass and fail directories not found; did basecalling fail?')

    fastq_fname = 'fastq_%s.tgz' % args.jobname
    if os.path.exists(fastq_fname):
      info('%s already exists; skipping' % fastq_fname)
    else:
      run(['tar', '-czvf', fastq_fname, 'pass', 'fail'])
      upload_public(fastq_fname)

  consensus_fname = 'consensus_%s.fasta' % args.jobname
  if os.path.exists(consensus_fname):
    info('consensus results already exist, skipping')
  else:
    run([
      os.path.join(THISDIR, 'callco.jl'),
      # size selection lower limit
      '--b1', '1000',
      # size selection upper limit
      '--b2', '7000',
      '--out', consensus_fname,
      'pass'])

    upload_public(consensus_fname)

  if args.upload_fastq:
    success('fastq results available at '
            'https://ontseq-res-public.s3.amazonaws.com/%s' % fastq_fname)

  success('consensus results available at '
          'https://ontseq-res-public.s3.amazonaws.com/%s' % consensus_fname)

if __name__ == '__main__':
  start()
