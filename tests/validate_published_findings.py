#!/usr/bin/env python3
import gzip,json,sys
errors=[];total=0
with gzip.open(sys.argv[1],'rt') as h:
 for total,line in enumerate(h,1):
  r=json.loads(line);a=int(r['ExactAltReads']);ref=int(r['CleanReferenceReads']);ex=int(r['ExcludedReads']);call=int(r.get('CallableReads',r['CallableAlignments']));exam=int(r.get('TotalReadsExamined',r['UniqueAlignments']))
  status='MIXED_ALT_AND_REFERENCE' if a and ref else 'ALT_SUPPORTED' if a else 'NO_EXACT_ALT_SUPPORT' if ref else 'NO_CALLABLE_READS'
  if call!=a+ref or exam!=call+ex or r['ReadValidationStatus']!=status:errors.append(r['EventID'])
  fraction=r.get('AltFractionAmongClean')
  if call==0 and fraction not in (None,'NA','N/A',''):errors.append(r['EventID']+':undefined_fraction')
if errors:raise SystemExit(f'FAIL: {len(errors)} invariant errors; first={errors[:10]}')
print(f'published finding validation: PASS ({total} records)')
