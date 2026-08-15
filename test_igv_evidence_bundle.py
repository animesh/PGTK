#!/usr/bin/env python3
import subprocess,sys,tempfile
from pathlib import Path
import pysam
root=Path(__file__).resolve().parent
with tempfile.TemporaryDirectory(prefix='pgtk_igv_') as td:
 t=Path(td);genome=t/'genome.fa';genome.write_text('>1\n'+'A'*500+'\n>2\n'+'C'*500+'\n>7\n'+'G'*500+'\n');pysam.faidx(str(genome))
 vcf=t/'TK1.rna.validated.vcf';vcf.write_text('##fileformat=VCFv4.2\n##contig=<ID=1,length=500>\n##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: Allele|SYMBOL|PICK|CANONICAL">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t100\t.\tA\tG\t.\tPASS\tCSQ=G|GENE1|1|YES\n')
 fusion=t/'TK1.fusion.tsv';fusion.write_text('#gene1\tgene2\tbreakpoint1\tbreakpoint2\nA\tB\t2:200\t7:300\n')
 splice=t/'TK1.splice.tsv';splice.write_text('Event\tStatus\tJunctions\nTX1\tRNA_VALIDATED\t1:150-180\n')
 bam=t/'TK1.Aligned.sortedByCoord.out.bam';header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'1','LN':500},{'SN':'2','LN':500},{'SN':'7','LN':500}]}
 with pysam.AlignmentFile(bam,'wb',header=header) as out:
  for rid,start,name in ((0,90,'r1'),(0,145,'r2'),(1,195,'r3'),(2,295,'r4')):
   r=pysam.AlignedSegment();r.query_name=name;r.reference_id=rid;r.reference_start=start;r.mapping_quality=60;r.cigarstring='20M';r.query_sequence='A'*20;r.query_qualities=pysam.qualitystring_to_array('I'*20);out.write(r)
 pysam.index(str(bam));prefix=t/'pgtk_igv'
 subprocess.run([sys.executable,str(root/'build_igv_evidence_bundle.py'),'--genome',str(genome),'--rna-vcf',str(vcf),'--fusion-table',str(fusion),'--splice-table',str(splice),'--bam',f'TK1={bam}','--output-prefix',str(prefix)],check=True)
 outbam=Path(str(prefix)+'.TK1.events.bam');assert outbam.is_file() and Path(str(outbam)+'.bai').is_file();assert pysam.quickcheck(str(outbam))==''
 with pysam.AlignmentFile(outbam,'rb') as h:assert h.count(until_eof=True)==4
 assert sum(1 for _ in open(str(prefix)+'.events.tsv'))==4
print('PASS: Pysam IGV bundle fixture')
