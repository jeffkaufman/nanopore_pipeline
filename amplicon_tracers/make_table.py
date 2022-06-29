#!/usr/bin/env python3

import sys
import json
from collections import defaultdict

from determine_matches import parse_sequencing_barcodes

def start(tracer_jsons, tracer_table):
    header = ['indexing_barcode_id']
    sequencing_barcodes = parse_sequencing_barcodes()
    for sequencing_barcode_id, _ in sequencing_barcodes:
        header.append(sequencing_barcode_id)
    output = [header]

    # idx barcode -> [sq 1 count, sq 2 count, ...]
    counts = defaultdict(lambda: [0]*len(sequencing_barcodes))

    with open(tracer_jsons) as inf:
        for line in inf:
            line = line.strip()
            if not line: continue

            sequencing_barcode, *scores = json.loads(line)
            # TODO: mark as unclassified if score is below a threshold
            counts[sequencing_barcode][scores.index(max(scores))] += 1

    for sequencing_barcode in sorted(counts):
        row = [sequencing_barcode]
        row.extend(counts[sequencing_barcode])
        output.append(row)

    with open(tracer_table, 'w') as outf:
        for row in output:
            outf.write('\t'.join([str(x) for x in row]))
            outf.write('\n')

if __name__ == "__main__":
    start(*sys.argv[1:])
