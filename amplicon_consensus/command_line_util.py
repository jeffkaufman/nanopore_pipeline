import subprocess
import os
import sys
import re

COLORS = {
  'red': '\033[91m',
  'green': '\033[92m',
  'magenta': '\033[95m',
  'cyan': '\033[96m',
  'end': '\033[0m',
}

THISDIR = os.path.abspath(os.path.dirname(__file__))

def print_color(color, msg, file=sys.stdout):
  print('%s%s%s' % (COLORS[color], msg, COLORS['end']), file=file)

def info(msg):
  print_color('cyan', msg)

def success(msg):
  print_color('green', msg)

def die(msg):
  print_color('red', msg)
  sys.exit(1)

def run(command, capture_output=False):
  print_color('magenta', '\n%s\n' % ' '.join(command))
  if capture_output:
    return subprocess.check_output(command)
  else:
    subprocess.check_call(command)

def validate_jobname(jobname):
  if not re.match(r'^[a-zA-Z0-9_-]+$', jobname):
    die('job names must be alphanumeric; got %r' % jobname)

def tarfile_name(jobname):
  return 'raw-%s.tar' % jobname
    
def s3_url(jobname):
  return 's3://summer-seq-data/%s' % tarfile_name(jobname)
