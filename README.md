# Nanopore Pipeline

## Amplicon Consensus

This pipeline takes sequencing results, calls bases, demuxes, and determines
consensus sequences.

### Pre-setup

1. Generate public key: `ssh-keygen -t ed25519 -C "your_email@example.com"`
1. Copy the public key: `cat ~/.ssh/ed25519.pub | pbcopy`
1. Add the key to github: https://github.com/settings/ssh/new

### Setup

1. Get Jeff or another existing user to add your public key to the EC2 instance
   (jefftk-basecaller) and to the sequencing desktop (summer)
   1. These should be appended to `~/.ssh/authorized_keys`
1. Verify you can ssh to the sequencing desktop (summer):
   `ssh summer@10.114.5.117 -p 64355`.  If you're on WiFi you need `MIT Secure`,
   not `MIT Guest` or `MIT`.  Not sure if this can be done when working from
   home.

### Operation

1. `nanopore_pipeline/amplicon_consensus/local.py <jobname> <bararr>`

Run `nanopore_pipeline/amplicon_consensus/local.py --help` and
`nanopore_pipeline/amplicon_consensus/call_bases_and_consense.py` for more
explanation and options.

### Troubleshooting

The program prints information about what it's running in blue and commands it
runs in magenta.

It untars the results to [jobname] (ex: 220603anj) and then works in
[jobname]Results (ex: 220603anjResults).

It avoids repeating work by checking if outputs already exist.  If it gets
that wrong, logging into the EC2 instance and deleting outputs will fix it.

File bugs at https://github.com/jeffkaufman/nanopore_pipeline/issues/new
