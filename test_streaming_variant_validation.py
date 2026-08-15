#!/usr/bin/env python3
import csv
import subprocess
import tempfile
from pathlib import Path
import pysam

root = Path(__file__).resolve().parent
with tempfile.TemporaryDirectory(prefix='pgtk_streaming_fixture_') as temporary:
    work = Path(temporary)
    fasta = work / 'genome.fa'
    fasta.write_text('>chr1\n' + 'A' * 500 + '\n', encoding='utf-8')
    pysam.faidx(str(fasta))
    header = {'HD': {'VN': '1.6', 'SO': 'coordinate'}, 'SQ': [{'SN': 'chr1', 'LN': 500}]}
    bam = work / 'TEST.bam'
    with pysam.AlignmentFile(bam, 'wb', header=header) as output:
        for index, base in enumerate(('T', 'T', 'T', 'A')):
            read = pysam.AlignedSegment()
            read.query_name = f'read{index}'
            read.query_sequence = 'A' * 49 + base + 'A' * 50
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 50
            read.mapping_quality = 60
            read.cigar = [(0, 100)]
            read.query_qualities = pysam.qualitystring_to_array('I' * 100)
            output.write(read)
    pysam.index(str(bam))
    vcf = work / 'TEST.rna.validated.vcf'
    vcf.write_text(
        '##fileformat=VCFv4.2\n'
        '##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: Allele|Consequence|SYMBOL|Feature|ENSP|HGVSc|HGVSp|Codons|Amino_acids">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        'chr1\t100\tv1\tA\tT\t.\tPASS\tCSQ=T|missense_variant|GENE|TX|P|c.1A>T|p.K*|AAA/TAA|K/*\n',
        encoding='utf-8',
    )
    vcfgz = Path(str(vcf) + '.gz')
    pysam.tabix_compress(str(vcf), str(vcfgz), force=True)
    pysam.tabix_index(str(vcfgz), preset='vcf', force=True)
    samples = work / 'samples.csv'
    samples.write_text('sample,srr\nTEST,SRRTEST\n', encoding='utf-8')
    codon_prefix = work / 'TEST.variant_codon_validation'
    subprocess.run([
        'python3', str(root / 'validate_variant_codons.py'), '--vcf', str(vcfgz),
        '--bam', f'TEST={bam}', '--genome', str(fasta), '--threads', '2',
        '--min-alt-reads', '3', '--min-alt-fraction', '0.05', '--output-prefix', str(codon_prefix),
    ], check=True)
    provenance_prefix = work / 'TEST.variant_read_provenance'
    subprocess.run([
        'python3', str(root / 'validate_variant_read_provenance.py'), '--vcf', str(vcfgz),
        '--bam', f'TEST={bam}', '--samples', str(samples), '--threads', '2',
        '--output-prefix', str(provenance_prefix),
    ], check=True)
    with Path(f'{codon_prefix}.all.tsv').open(encoding='utf-8') as handle:
        codon_rows = list(csv.DictReader(handle, delimiter='\t'))
    with Path(f'{provenance_prefix}.supporting_reads.tsv').open(encoding='utf-8') as handle:
        provenance_rows = list(csv.DictReader(handle, delimiter='\t'))
    assert len(codon_rows) == 1
    assert codon_rows[0]['ALT-supporting reads'] == '3'
    assert codon_rows[0]['Overall status'] == 'VARIANT_CODON_VALIDATED'
    assert len(provenance_rows) == 3
    assert all(row['Observed allele'] == 'T' for row in provenance_rows)
print('PASS: indexed streaming codon and provenance validation fixture')
