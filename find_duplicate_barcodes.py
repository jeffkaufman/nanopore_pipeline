import sys
from collections import defaultdict
import argparse

def start():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'bararr', help='filename for barcode metadata in tab-delimited format.')
    args = parser.parse_args()

    locs = defaultdict(list)

    with open(args.bararr) as inf:
        for index, line in enumerate(inf):
            if index == 0: continue  # skip headers

            line = line.strip()
            if not line: continue

            sample_id, _, fwd, rev = line.split('\t')

            locs[fwd].append('%s:%s:fwd' % (sample_id, index))
            locs[rev].append('%s:%s:rev' % (sample_id, index))

    for seq, seq_locs in locs.items():
        if len(seq_locs) != 1:
            print('%s: %s' % (seq, ', '.join(seq_locs)))

if __name__ == '__main__':
    start()
