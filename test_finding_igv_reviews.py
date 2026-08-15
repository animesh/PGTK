#!/usr/bin/env python3
import importlib.util
from pathlib import Path

root = Path(__file__).resolve().parent

def load(name, filename):
    spec = importlib.util.spec_from_file_location(name, root / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

reviews = load('finding_reviews', 'build_finding_igv_reviews.py')
bundle = load('igv_bundle', 'build_igv_evidence_bundle.py')

def alignment(sequence, cigar, start=100, quality='I'):
    return ['read1', '0', '1', str(start), '60', cigar, '*', '0', '0', sequence, quality * len(sequence)]

# SNP, normalized-anchor insertion, adjacent repeat-equivalent insertion, deletion.
assert reviews.classify(alignment('G', '1M'), 100, 'A', 'G', 95, 110, 20, 20)[0] == 'EXACT_ALT'
assert reviews.classify(alignment('ACA', '1M2I'), 100, 'A', 'ACA', 95, 110, 20, 20)[0] == 'EXACT_ALT'
assert reviews.classify(alignment('ACA', '2M1I'), 100, 'A', 'AA', 95, 110, 20, 20)[0] == 'EXACT_ALT'
assert reviews.classify(alignment('AC', '1M1D1M'), 100, 'AT', 'A', 95, 110, 20, 20)[0] == 'EXACT_ALT'

# Gene identifiers must be searchable and sample-specific.
event = {'Label':'C>CA','Gene':'DDX1','Sample':'TK14','Chrom':'2','Start0':'15597433'}
assert reviews.safe_id(event).startswith('DDX1_TK14_2_15597434_C_CA')

# VEP selection must be ALT-specific and prefer PICK over canonical/first-record order.
annotations = [
    {'Allele':'G','SYMBOL':'WRONG','CANONICAL':'YES'},
    {'Allele':'A','SYMBOL':'DDX1','PICK':'1'},
]
assert bundle.choose_annotation(annotations, 'A', 'C')['SYMBOL'] == 'DDX1'
assert bundle.choose_annotation(annotations, 'G', 'C')['SYMBOL'] == 'WRONG'

# Every generated IGV batch must terminate after snapshot generation.
source = (root / 'build_finding_igv_reviews.py').read_text()
assert "snapshot {event_id}.png\\nexit\\n" in source
print('PASS: finding IGV classifier, gene naming, ALT-specific VEP selection and batch termination')
