#!/usr/bin/env python3
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from variant_read_evidence import classify_sam_fields,evidence_status
assert evidence_status(0,0,9)[3] is None
assert evidence_status(2,0,0)[0]=='ALT_SUPPORTED'
assert evidence_status(2,3,0)[0]=='MIXED_ALT_AND_REFERENCE'
base=['r','0','1','100','60','10M','*','0','0','ACGTACGTAC','IIIIIIIIII']
assert classify_sam_fields(base,101,'C','T',20,20)[0]=='CLEAN_REFERENCE'
mnv=base.copy();mnv[9]='ATGTACGTAC';assert classify_sam_fields(mnv,101,'CG','TG',20,20)[0]=='EXACT_ALT'
ins=['r','0','1','100','60','2M2I8M','*','0','0','ACGGGTACGTAC','IIIIIIIIIIII'];assert classify_sam_fields(ins,101,'C','CGG',20,20)[0]=='EXACT_ALT'
delr=['r','0','1','100','60','2M2D8M','*','0','0','ACGTACGTAC','IIIIIIIIII'];assert classify_sam_fields(delr,101,'CGT','C',20,20)[0]=='EXACT_ALT'
print('semantic regression tests: PASS')
