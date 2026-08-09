#!/usr/bin/env python3
import csv, subprocess, sys, tempfile
from pathlib import Path

SCRIPT=Path(__file__).with_name('expression_go_analysis.py')

def run(*args): subprocess.run([sys.executable,str(SCRIPT),*map(str,args)],check=True)
def rows(path):
    with open(path,encoding='utf-8',newline='') as h: return list(csv.DictReader(h,delimiter='\t'))

with tempfile.TemporaryDirectory() as tmp:
    d=Path(tmp); gtf=d/'genes.gtf'; go=d/'go.tsv'; samples=d/'samples.csv'
    with gtf.open('w') as h:
        for i in range(1,11): h.write(f'chr1\ttest\texon\t{i*100}\t{i*100+99}\t.\t+\t.\tgene_id "g{i}"; gene_name "G{i}"; gene_biotype "protein_coding";\n')
    go.write_text('Gene\tGO_ID\tGO_Name\tNamespace\n' + ''.join(f'G{i}\tGO:1\tPathway one\tbiological_process\n' for i in range(1,6)) + ''.join(f'G{i}\tGO:2\tPathway two\tbiological_process\n' for i in range(6,11)))
    samples.write_text('TK,Group,sample,fastq_1,fastq_2,baseline,CDS\nP1,base,S1,a,b,true,c\nP1,prog,S2,d,e,false,f\n')
    for sample,vals in [('S1',[10,10,10,10,10,1,1,1,1,1]),('S2',[1,1,1,1,1,10,10,10,10,10])]:
        out=d/f'{sample}.gene_counts.tsv'
        out.write_text('# Program:featureCounts\nGeneid\tChr\tStart\tEnd\tStrand\tLength\tbam\n'+''.join(f'g{i}\tchr1\t{i*100}\t{i*100+99}\t+\t100\t{v}\n' for i,v in enumerate(vals,1)))
    run('merge-counts','--counts',d/'S1.gene_counts.tsv',d/'S2.gene_counts.tsv','--gtf',gtf,'--output-prefix',d/'merged')
    matrix=d/'merged.gene_expression.tsv'; assert len(rows(matrix))==10
    assert abs(sum(float(r['S1_TPM']) for r in rows(matrix))-1e6)<1e-6
    run('expression-go','--matrix',matrix,'--samples',samples,'--go-mapping',go,'--go-min-size','2','--go-max-size','20','--cpm-threshold','1','--output-prefix',d/'ego')
    assert rows(d/'ego.expression_ora.tsv') and rows(d/'ego.ranked_go.tsv')
    assert all(r['FDRThreshold']=='0.1' for r in rows(d/'ego.summary.tsv'))
    run('sample-ora','--matrix',matrix,'--sample','S1','--subject','P1','--group','base','--all-samples','S1,S2','--go-mapping',go,'--go-min-size','2','--go-max-size','20','--output-prefix',d/'S1_ora')
    run('ranked-go','--matrix',matrix,'--sample','S2','--baseline-sample','S1','--subject','P1','--group','prog','--all-samples','S1,S2','--go-mapping',go,'--go-min-size','2','--go-max-size','20','--output-prefix',d/'S2_ranked')
    ranked_rows=rows(d/'S2_ranked.ranked_go.tsv')
    ranked_by_go={r['GO_ID']:r for r in ranked_rows}
    assert float(ranked_by_go['GO:1']['MeanScore']) < 0
    assert float(ranked_by_go['GO:2']['MeanScore']) > 0
    ranked_summary=rows(d/'S2_ranked.summary.tsv')[0]
    assert ranked_summary['Analysis']=='S2_vs_S1:log2_tpm_fold_change'
    assert ranked_summary['Sample']=='S2' and ranked_summary['BaselineSample']=='S1'
    assert int(ranked_summary['NonZeroScores'])==10
    assert int(ranked_summary['PositiveScores'])==5
    assert int(ranked_summary['NegativeScores'])==5
    assert float(ranked_summary['MinScore']) < 0 < float(ranked_summary['MaxScore'])
    assert ranked_summary['Pseudocount']=='0.5'
    run('merge-expression-go','--samples',samples,'--ora',d/'S1_ora.expression_ora.tsv','--ranked',d/'S2_ranked.ranked_go.tsv','--summary',d/'S1_ora.summary.tsv',d/'S2_ranked.summary.tsv','--output-prefix',d/'scatter')
    assert rows(d/'scatter.expression_ora.tsv') and rows(d/'scatter.ranked_go.tsv')
    self_compare=subprocess.run([sys.executable,str(SCRIPT),'ranked-go','--matrix',str(matrix),'--sample','S1','--baseline-sample','S1','--subject','P1','--group','base','--all-samples','S1,S2','--go-mapping',str(go),'--go-min-size','2','--go-max-size','20','--output-prefix',str(d/'self')],text=True,capture_output=True)
    assert self_compare.returncode != 0 and '--sample and --baseline-sample must be different' in self_compare.stderr
    missing_compare=subprocess.run([sys.executable,str(SCRIPT),'ranked-go','--matrix',str(matrix),'--sample','MISSING','--baseline-sample','S1','--subject','P1','--group','prog','--all-samples','S1,S2','--go-mapping',str(go),'--go-min-size','2','--go-max-size','20','--output-prefix',str(d/'missing')],text=True,capture_output=True)
    assert missing_compare.returncode != 0 and 'MISSING_TPM' in missing_compare.stderr
    identical=d/'identical.tsv'
    identical.write_text('Gene\tA_raw_count\tA_TPM\tB_raw_count\tB_TPM\nG1\t1\t2\t1\t2\nG2\t1\t3\t1\t3\n')
    zero_compare=subprocess.run([sys.executable,str(SCRIPT),'ranked-go','--matrix',str(identical),'--sample','B','--baseline-sample','A','--subject','P1','--group','prog','--all-samples','A,B','--go-mapping',str(go),'--go-min-size','1','--go-max-size','20','--output-prefix',str(d/'zero')],text=True,capture_output=True)
    assert zero_compare.returncode != 0 and '0 non-zero scores' in zero_compare.stderr
    invalid=d/'invalid.csv'; invalid.write_text('sample,TK,Group,baseline\nA,P2,x,false\nB,P2,y,false\n')
    run('merge-expression-go','--samples',invalid,'--output-prefix',d/'invalid_merge')
    skipped=rows(d/'invalid_merge.summary.tsv')
    assert len(skipped)==2 and all(r['Status']=='SKIPPED' and 'no baseline' in r['Message'] for r in skipped)
    variant_samples=d/'variant_samples.csv'
    variant_samples.write_text('sample,TK,Group,baseline\nS1,P1,base,true\nS2,P1,prog_a,false\nS3,P1,prog_b,false\n')
    for sample,genes in [('S2',['G1','G2','G7','G8']),('S3',['G1','G3','G8','G9'])]:
        (d/f'{sample}.genes.tsv').write_text('Sample\tGene\n'+''.join(f'{sample}\t{g}\n' for g in genes))
    run('variant-sets','--genes',d/'S2.genes.tsv',d/'S3.genes.tsv','--samples',variant_samples,'--go-mapping',go,'--gtf',gtf,'--biotypes','protein_coding','--go-min-size','2','--go-max-size','20','--output-prefix',d/'vgo')
    summary=rows(d/'vgo.summary.tsv'); labels={r['Analysis'] for r in summary}
    assert labels == {'common_all_progression','S2_all','S2_exclusive','S3_all','S3_exclusive'}
    assert all(r['FDRThreshold']=='0.1' for r in summary)
    by_label={r['Analysis']:r for r in summary}
    assert by_label['S2_all']['ForegroundGenes']=='4'
    assert by_label['S3_all']['ForegroundGenes']=='4'
    assert by_label['common_all_progression']['ForegroundGenes']=='2'
    run('variant-sets','--genes',d/'S2.genes.tsv',d/'S3.genes.tsv','--samples',variant_samples,'--go-mapping',go,'--gtf',gtf,'--biotypes','protein_coding','--go-min-size','2','--go-max-size','20','--fdr-threshold','0','--output-prefix',d/'vgo_zero')
    assert all(r['FDRThreshold']=='0.0' and r['SignificantGOTerms']=='0' for r in rows(d/'vgo_zero.summary.tsv'))
print('expression and progression-set GO tests passed')
