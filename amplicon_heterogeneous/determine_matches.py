import glob
from collections import namedtuple
from collections import defaultdict
from Bio import Align
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO import SeqRecord
import json
import sys

Barcode = namedtuple('Barcode', ['barcode_id', 'fwd', 'rev'])

def parse_barcodes():
    barcodes = {}
    with open("bararr.csv") as inf:
        for i, line in enumerate(inf):
            if i == 0:
                continue
            line = line.strip()
            cols = line.split("\t")
            if len(cols) == 4:
                barcode_id, _, fwd, rev = cols
                guppy_barcode_id = 'bw_%s' % (str(len(barcodes)).zfill(3))
                barcodes[guppy_barcode_id] = Barcode(
                    barcode_id, Seq(fwd), Seq(rev))

    return barcodes

def parse_tracers():
    return [
        ["tr1", "AGAGGCTCTAAGGGCTTCTCAGTGTATGTTACATCCC"],
        ["tr2", "AGCGCCTAACGCCGCTCCTTTCGCTTTCTTCCCTTT"],
    ] # FIXME these are made up

def start(fnames):
    barcodes = parse_barcodes()
    tracers = parse_tracers()

    aligner = Align.PairwiseAligner(
        # These are the scoring settings porechop uses by default.
        # https://github.com/rrwick/Porechop/blob/master/porechop/porechop.py#L145
        end_gap_score = 0,
        match_score = 3,
        mismatch_score = -6,
        internal_open_gap_score = -5,
        internal_extend_gap_score = -2)

    for fname in fnames:
        with open(fname) as inf:
            records = SeqIO.parse(inf, 'fastq')
            for record in records:
                handle_record(fname, record, barcodes, tracers, aligner)

def orient(record, barcode, aligner):
    fwd_score = aligner.score(record.seq, barcode.fwd) + aligner.score(
        record.seq, barcode.rev.reverse_complement())
    rev_score = aligner.score(record.seq, barcode.rev) + aligner.score(
        record.seq, barcode.fwd.reverse_complement())

    if rev_score > fwd_score:
        record.seq = record.seq.reverse_complement()

def trim_barcodes(record, barcode, aligner):
    seq = record.seq

    fwd_alignment = aligner.align(seq, barcode.fwd)[0]
    # Look at the aligned spans, take the last one, and then take the end of
    # that span.
    seq = record.seq[fwd_alignment.aligned[0][-1][-1]:]
    rev_alignment = aligner.align(seq, barcode.rev.reverse_complement())[0]
    # Look at the aligned spans, take the first one, and then take the
    # beginning of that span.
    seq = seq[:rev_alignment.aligned[0][0][0]]
    return seq

def handle_record(fname, record, barcodes, tracers, aligner):
    guppy_barcode = record.description.split()[-1].removeprefix("barcode=")
    barcode = barcodes[guppy_barcode]

    orient(record, barcode, aligner)
    seq = trim_barcodes(record, barcode, aligner)
    if not seq:
        return

    result = [barcode.barcode_id]
    for tracer_id, tracer_seq in tracers:
        result.append(int(aligner.score(seq, tracer_seq)))
    print(json.dumps(result))

if __name__ == "__main__":
    start(sys.argv[1:])
