#!/usr/bin/python3

import sys
import argparse
import re
import time
import os

from command_line_util import *

def run_on_sequencing_desktop(cmd, **kwargs):
  return run(['ssh', 'summer@10.114.0.70', '-p', '64355', cmd], **kwargs)

def s3_contains_object(s3_url):
  try:
    run_on_sequencing_desktop('aws s3 ls %s > /dev/null' % s3_url)
    return True
  except subprocess.CalledProcessError:
    return False

def copy_raw_data_to_s3(jobname, force_upload):
  if s3_contains_object(s3_url(jobname)):
    if force_upload:
      info('%s already on s3, but overwriting because of --force-upload' %
           jobname)
    else:
      info('%s already on s3, skipping' % jobname)
      return

  run_on_sequencing_desktop(
    'set -o pipefail && cd minknow/data/ && tar -cf - %s |'
    ' aws s3 cp - %s' % (jobname, s3_url(jobname)))

def parse_aws_output(output, key, fallback=None):
  output = (output or b'').decode('utf-8')

  for line in output.split('\n'):
    if ' %s ' % key in line:
      expected = ''
      try:
        _, expected, _, value, _ = line.split()
      except ValueError:
        pass
      if expected != key:
        print(output)
        die('unexpected output: %r' % line)
      return value
  if fallback:
    return fallback

  print(output)
  die('unable to parse aws output')

# Possible known outputs:
#   stopped
#   initializing
#   ok
def ec2_status():
  output = run_on_sequencing_desktop(
    'aws ec2 describe-instance-status --instance-id i-0213be8d87b18ef2c',
    capture_output=True)
  # print(output.decode('utf-8'))
  return parse_aws_output(output, 'Status', fallback='stopped')

def start_ec2():
  status = ec2_status()
  if status == 'stopped':
    info('starting ec2 instance...')
    run_on_sequencing_desktop(
      'aws ec2 start-instances --instance-id i-0213be8d87b18ef2c')
    time.sleep(60)
    status = ec2_status()

  while status == 'initializing':
    info('waiting for ec2 instance to finish initializing...')
    time.sleep(5)
    status = ec2_status()

  if status != 'ok':
    die('failed to start ec2 instance')

  info('ec2 instance is up')

def ec2_public_dns():
  output = run_on_sequencing_desktop(
    'aws ec2 describe-instances --instance-id i-0213be8d87b18ef2c',
    capture_output=True)
  return parse_aws_output(output, 'PublicDnsName')

def ec2_address():
  return 'ec2-user@%s' % ec2_public_dns()

def run_on_ec2(cmd, **kwargs):
  return run(['ssh', ec2_address(), cmd], **kwargs)

def process_on_ec2(args):
  remote_indexing_barcodes_name = '%s.bararr' % args.jobname
  run(['scp', args.indexing_barcodes, '%s:%s' % (
    ec2_address(), remote_indexing_barcodes_name)])

  remote_sequencing_barcodes_name = None
  if args.sequencing_barcodes:
    remote_sequencing_barcodes_name = '%s.sequencing_barcodes' % args.jobname
    run(['scp', args.sequencing_barcodes, '%s:%s' % (
      ec2_address(), remote_sequencing_barcodes_name)])
    args.remote_args += (
      ' --sequencing-barcodes %s' % remote_sequencing_barcodes_name)

  run_on_ec2(
    'python3 nanopore_pipeline/remote.py'
    ' %s %s %s' % (args.jobname, remote_indexing_barcodes_name, args.remote_args))

def stop_ec2():
  info("stopping ec2 instance...")
  output = run_on_sequencing_desktop(
    'aws ec2 stop-instances --instance-id i-0213be8d87b18ef2c',
    capture_output=True)
  if parse_aws_output(output, '64'):
    success('instance is stopping')

def start():
  parser = argparse.ArgumentParser(
    description='Automate sequencing end to end')
  parser.add_argument(
    '--remote-args',
    help='Arguments to pass through to the the ec2 instance.  Run '
    '"remote.py --help" to see available options.  Ex: '
    '--remote-args="--upload-fastq --max-length=10000"')
  parser.add_argument(
    '--leave-ec2-running', action='store_true',
    help='By default this script shuts down the ec2 instance when done, but '
    'if you are iterating you probably want to leave it up until you are '
    'finished.')
  parser.add_argument(
    '--force-upload', action='store_true',
    help='Upload to s3 even if there are already results uploaded under this '
    'name.  Use when incrementally processing sequencer output.')
  parser.add_argument(
    '--sequencing-barcodes', metavar='FNAME',
    help='If present, use these barcodes to determine distribution of tracers')
  parser.add_argument('--determine-consensus', action='store_true')
  parser.add_argument('jobname', help='experiment name, ex: 220603anj')
  parser.add_argument(
    'indexing_barcodes',
    help='filename for barcode metadata in tab-delimited format.')
  args = parser.parse_args()

  validate_jobname(args.jobname)

  if not os.path.exists(args.indexing_barcodes):
    die('expected %r to exist' % args.indexing_barcodes)

  args.remote_args = args.remote_args or ''

  if args.force_upload:
    args.remote_args += ' --force-download'
  if args.determine_consensus:
    args.remote_args += ' --determine-consensus'

  copy_raw_data_to_s3(args.jobname, args.force_upload)
  try:
    start_ec2()
    process_on_ec2(args)
  finally:
    if not args.leave_ec2_running:
      stop_ec2()

if __name__ == '__main__':
  start()
