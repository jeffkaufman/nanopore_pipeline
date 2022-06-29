
#!/usr/bin/python3

import sys
import json
import glob
import os
import subprocess
import shutil
import argparse
import filecmp

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
    '-c', os.path.join(os.path.expanduser('~/ont-guppy/data'), basecall_model),
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

def validate_consensus(consensus_fname):
  n_bad = 0
  n_total = 0
  with open(consensus_fname) as inf:
    for line in inf:
      if line.startswith('>'):
        n_total += 1
      elif line.strip() in ['AAAA', 'CCCC']:
        n_bad += 1
  if n_bad == n_total:
    die('Consensus failed.  Do you need to tweak --min-length or '
        '--max-length?')
  info('Consensus found results for %s/%s' % (n_total - n_bad, n_total))

def start():
  parser = argparse.ArgumentParser(
    description='Call bases and determine consensus sequences')
  parser.add_argument(
    '--upload-fastq', action='store_true',
    help='Upload FASTQ files to S3 after basecalling.')
  parser.add_argument(
    '--min-length', type=int,
    default=1000,
    help='Minimum sequence length to identify during consensus')
  parser.add_argument(
    '--max-length', type=int,
    default=7000,
    help='Maximum sequence length to identify during consensus')
  parser.add_argument(
    '--force-download', action='store_true',
    help='Download from s3 even if there is already a tarball present.')
  parser.add_argument('--sequencing-barcodes', metavar='FNAME')
  parser.add_argument('--determine-consensus', action='store_true')

  parser.add_argument('jobname')
  parser.add_argument('indexing_barcodes')
  args = parser.parse_args()

  validate_jobname(args.jobname)
  if not os.path.exists(args.indexing_barcodes):
    die('expected %r to exist' % args.indexing_barcodes)

  indexing_barcodes_in = os.path.abspath(args.indexing_barcodes)

  tarfile = tarfile_name(args.jobname)

  tarfile_prev = tarfile + '.prev'
  if not os.path.exists(tarfile) or args.force_download:
    if os.path.exists(tarfile):
      os.replace(tarfile, tarfile_prev)

    run(['aws', 's3', 'cp',  s3_url(args.jobname), '.'])

  raw_dir = os.path.abspath(args.jobname)
  results_dir = os.path.abspath('%sResults' % args.jobname)

  if os.path.exists(tarfile_prev):
    if not filecmp.cmp(tarfile, tarfile_prev):
      info('tarfile changed, will redo any steps from before')
      if os.path.exists(raw_dir):
        shutil.rmtree(raw_dir)
      if os.path.exists(results_dir):
        shutil.rmtree(results_dir)

  if not os.path.exists(raw_dir):
    os.mkdir(raw_dir)
    run(['tar', '-xvf', tarfile, '-C', raw_dir])

  indexing_barcodes_in_results = os.path.join(results_dir, 'bararr.csv')
  if os.path.exists(results_dir) and (
      not os.path.exists(indexing_barcodes_in_results) or not filecmp.cmp(
        indexing_barcodes_in, indexing_barcodes_in_results)):
    info('new indexing_barcodes')
    shutil.rmtree(results_dir)

  sequencing_barcodes_in_results = os.path.join(
    results_dir, 'sequencing_barcodes')
  tracer_table = os.path.join(results_dir, 'tracer_table.tsv')
  if args.sequencing_barcodes:
    if not os.path.exists(sequencing_barcodes_in_results) or not filecmp.cmp(
        args.sequencing_barcodes, sequencing_barcodes_in_results):
      info('new tracer ids')
      shutil.copyfile(args.sequencing_barcodes, sequencing_barcodes_in_results)
      if os.path.exists(tracer_table):
        os.remove(tracer_table)

  if not os.path.exists(results_dir):
    os.mkdir(results_dir)
    shutil.copyfile(indexing_barcodes_in, indexing_barcodes_in_results)

  os.chdir(results_dir)
  run([
    os.path.join(THISDIR, 'bc_csv_to_fasta.py'),
    results_dir,
  ])  # creates <results_dir>/samples.csv

  if not os.path.exists('pass'):
    try:
      report_fname, = glob.glob(os.path.join(raw_dir, '**/report_*.json'),
                                recursive=True)
    except ValueError:
      die('No report*.json file -- has this run completed?')

    run_dir = os.path.dirname(report_fname)

    call_bases(run_dir, results_dir)

  result_dirs = [x for x in glob.glob('pass/*')]
  if not result_dirs or result_dirs == ['pass/unclassified']:
    # so we don't think basecalling succeeded
    os.remove(indexing_barcodes_in_results)
    die('demuxing failed: all outputs are unclassified.  Is this the right '
        'barcode metadata file?')

  if args.upload_fastq:
    upload_fastq()

  if args.determine_consensus:
    determine_consensus(args)

  if args.sequencing_barcodes:
    determine_tracer_composition(results_dir, tracer_table)

def upload_fastq():
    if not os.path.exists('pass') and not os.path.exists('fail'):
      die('pass and fail directories not found; did basecalling fail?')

    fastq_fname = 'fastq_%s.tgz' % args.jobname
    if os.path.exists(fastq_fname):
      info('%s already exists; skipping' % fastq_fname)
    else:
      run(['tar', '-czvf', fastq_fname, 'pass', 'fail'])
      upload_public(fastq_fname)

    success('fastq results available at '
            'https://ontseq-res-public.s3.amazonaws.com/%s' % fastq_fname)

def determine_consensus(args):
  consensus_fname = 'consensus_%s.fasta' % args.jobname

  min_length_fname = 'min_length'
  max_length_fname = 'max_length'
  if os.path.exists(consensus_fname) and os.path.exists(min_length_fname):
      with open(min_length_fname) as inf:
        if inf.read().strip() != str(args.min_length):
          info('min_length changed, discarding previous consensus')
          os.remove(consensus_fname)

  if os.path.exists(consensus_fname) and os.path.exists(max_length_fname):
    if os.path.exists(max_length_fname):
      with open(max_length_fname) as inf:
        if inf.read().strip() != str(args.max_length):
          info('max_length changed, discarding previous consensus')
          os.remove(consensus_fname)

  if os.path.exists(consensus_fname):
    info('consensus results already exist, skipping')
  else:
    with open(min_length_fname, 'w') as outf:
      outf.write(str(args.min_length))
    with open(max_length_fname, 'w') as outf:
      outf.write(str(args.max_length))

    run([
      os.path.join(THISDIR, 'amplicon_consensus', 'callco.jl'),
      # size selection lower limit
      '--b1', str(args.min_length),
      # size selection upper limit
      '--b2', str(args.max_length),
      '--out', consensus_fname,
      'pass'])

    validate_consensus(consensus_fname)

    upload_public(consensus_fname)

  success('consensus results available at '
          'https://ontseq-res-public.s3.amazonaws.com/%s' % consensus_fname)

def determine_tracer_composition(results_dir, tracer_table):
  tracer_jsons = os.path.join(results_dir, 'tracer_scores.jsons')
  if not os.path.exists(tracer_table):
    run([
      os.path.join(THISDIR, 'amplicon_tracers/determine_matches.py'),
      tracer_jsons])
    run([
      os.path.join(THISDIR, 'amplicon_tracers/make_table.py'),
      tracer_jsons, tracer_table])
    # TODO put on s3
  success('tracers processed')

if __name__ == '__main__':
  start()
