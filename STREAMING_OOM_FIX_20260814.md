# Streaming OOM fix 2026-08-14

This release is based on the exact Saga files captured from failed job 19195607, merged into the complete clean pipeline archive.

Changes:

- Replaced `pysam.samtools.view(..., catch_stdout=True)` in codon validation with indexed, per-variant `pysam.AlignmentFile.fetch()` streaming.
- Replaced the same whole-SAM capture pattern in read-provenance validation with indexed streaming.
- Preserved duplicate, secondary, supplementary, unmapped, MAPQ, base-quality, CIGAR and allele filters.
- Added BAM contig resolution for `chr` and non-`chr` naming.
- Set initial and retry memory for codon validation to 8, 16 and 32 GB.
- Set initial and retry memory for read provenance to 8, 16 and 32 GB.
- Set PyPGATK FASTA initial and retry memory to 8, 16 and 32 GB.
- Preserved the clean no-global-bind runtime, container-native Python, Pysam BGZF/tabix handling, and offline IGV Reports fixes.
- Added `test_streaming_variant_validation.py` to test codon counts and read provenance from a synthetic indexed BAM.
- Updated the pipeline validator to reject the old in-memory SAM-capture architecture.

The pipeline resumes from the existing Nextflow work directory. Successful upstream tasks remain cacheable; changed codon, provenance and dependent tasks rerun.
