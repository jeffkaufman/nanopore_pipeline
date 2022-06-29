#!/usr/bin/env python3
import glob
from collections import namedtuple
from collections import defaultdict
from Bio import Align
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO import SeqRecord
import json
import sys

IndexingBarcode = namedtuple(
    'IndexingBarcode', ['indexing_barcode_id', 'fwd', 'rev'])

def parse_indexing_barcodes():
    indexing_barcodes = {}
    with open("bararr.csv") as inf:
        for i, line in enumerate(inf):
            if i == 0: continue
            line = line.strip()
            if not line: continue

            indexing_barcode_id, _, fwd, rev = line.split("\t")
            guppy_indexing_barcode_id = 'bw_%s' % (
                str(len(indexing_barcodes)).zfill(3))
            indexing_barcodes[guppy_indexing_barcode_id] = IndexingBarcode(
                indexing_barcode_id, Seq(fwd), Seq(rev))

    return indexing_barcodes

def parse_sequencing_barcodes():
    sequencing_barcodes = []
    with open('sequencing_barcodes') as inf:
        for line in inf:
            line = line.strip()
            if line:
                sequencing_barcodes.append(line.split())
    return sequencing_barcodes

def start(tracer_jsons_fname):
    indexing_barcodes = parse_indexing_barcodes()
    sequencing_barcodes = parse_sequencing_barcodes()

    aligner = Align.PairwiseAligner(
        # These are the scoring settings porechop uses by default.
        # https://github.com/rrwick/Porechop/blob/master/porechop/porechop.py#L145
        end_gap_score = 0,
        match_score = 3,
        mismatch_score = -6,
        internal_open_gap_score = -5,
        internal_extend_gap_score = -2)

    with open(tracer_jsons_fname, 'w') as outf:
        for fname in glob.glob('pass/bw*/*.fastq'):
            with open(fname) as inf:
                records = SeqIO.parse(inf, 'fastq')
                for record in records:
                    handle_record(
                        fname, record, indexing_barcodes, sequencing_barcodes,
                        aligner, outf)

def orient(record, indexing_barcode, aligner):
    fwd_score = aligner.score(record.seq, indexing_barcode.fwd) + aligner.score(
        record.seq, indexing_barcode.rev.reverse_complement())
    rev_score = aligner.score(record.seq, indexing_barcode.rev) + aligner.score(
        record.seq, indexing_barcode.fwd.reverse_complement())

    if rev_score > fwd_score:
        record.seq = record.seq.reverse_complement()

def trim_indexing_barcodes(record, indexing_barcode, aligner):
    seq = record.seq

    fwd_alignment = aligner.align(seq, indexing_barcode.fwd)[0]
    # Look at the aligned spans, take the last one, and then take the end of
    # that span.
    seq = record.seq[fwd_alignment.aligned[0][-1][-1]:]
    rev_alignment = aligner.align(
        seq, indexing_barcode.rev.reverse_complement())[0]
    # Look at the aligned spans, take the first one, and then take the
    # beginning of that span.
    seq = seq[:rev_alignment.aligned[0][0][0]]
    return seq

def handle_record(fname, record, indexing_barcodes, sequencing_barcodes,
                  aligner, outf):
    guppy_indexing_barcode = record.description.split()[-1]
    if not guppy_indexing_barcode.startswith("barcode="):
        raise Exception("invalid guppy indexing barcode tag %r" %
                        guppy_indexing_barcode)
    guppy_indexing_barcode = guppy_indexing_barcode[len("barcode="):]
    indexing_barcode = indexing_barcodes[guppy_indexing_barcode]

    orient(record, indexing_barcode, aligner)
    seq = trim_indexing_barcodes(record, indexing_barcode, aligner)
    if not seq:
        return

    result = [indexing_barcode.indexing_barcode_id]
    for _, sequencing_barcode_seq in sequencing_barcodes:
        result.append(int(aligner.score(seq, sequencing_barcode_seq)))
    outf.write('%s\n' % json.dumps(result))

if __name__ == "__main__":
    start(*sys.argv[1:])
