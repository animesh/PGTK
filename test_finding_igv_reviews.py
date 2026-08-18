#!/usr/bin/env python3
import csv,importlib.util,subprocess,sys,tempfile
from pathlib import Path
root=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('reviews',root/'build_finding_igv_reviews.py'); reviews=importlib.util.module_from_spec(spec); spec.loader.exec_module(reviews)
def alignment(sequence,cigar,start=100):return ['read','0','1',str(start),'60',cigar,'*','0','0',sequence,'I'*len(sequence)]
assert reviews.classify(alignment('G','1M'),100,'A','G',95,110,20,20)[0]=='EXACT_ALT'
assert reviews.classify(alignment('ACA','1M2I'),100,'A','ACA',95,110,20,20)[0]=='EXACT_ALT'
assert reviews.classify(alignment('AC','1M1D1M'),100,'AT','A',95,110,20,20)[0]=='EXACT_ALT'
assert reviews.safe_id({'Label':'C>CA','Gene':'DDX1','Sample':'TK14','Chrom':'2','Start0':'15597433'}).startswith('DDX1_TK14_2_15597434_C_CA')
with tempfile.TemporaryDirectory() as directory:
    path=Path(directory); events=path/'events.tsv'; fields=['Event','Sample','Class','Chrom','Start0','End','Chrom2','Start2_0','End2','Label','Source','Gene','Consequence','Impact','Transcript','ProteinChange']; rows=[['E1','TK1','rna_variant','1','99','100','','','','A>G','a','G','','','',''],['E2','TK1','progression_variant','1','99','100','','','','A>G','b','G','','','',''],['E3','TK1','fusion','2','10','11','3','20','21','A--B','c','','','','','']]
    with events.open('w',newline='') as handle:writer=csv.writer(handle,delimiter='\t');writer.writerow(fields);writer.writerows(rows)
    output=path/'out';subprocess.run([sys.executable,str(root/'build_finding_igv_reviews.py'),'--events',str(events),'--genome','genome.fa','--output-dir',str(output),'--plan-only'],check=True)
    summary=(output/'consolidation_summary.txt').read_text();assert 'Consolidated findings: 2' in summary and 'Merged duplicate rows: 1' in summary;assert not any(item.is_dir() for item in output.iterdir())
source=(root/'build_finding_igv_reviews.py').read_text();assert 'BAM_CATEGORIES' in source and 'priority_findings.bed' in source
print('PASS: consolidated review classifier, deduplication, and flat output layout')
