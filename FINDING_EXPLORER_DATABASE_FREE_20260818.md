# Database-free complete finding explorer

The pipeline contains no alignment database and no metadata database. The failed monolithic IGV Reports stage and its finding-limit parameter were removed.

Every finding is written exactly once to a gzip-compressed JSON Lines partition determined by sample, evidence class, and chromosome. The local explorer loads these metadata partitions into memory, provides search and facets, and generates a standalone IGV.js report only for the selected locus. Generated locus reports are cached as ordinary HTML files.

Pipeline validation enforces:

- full manifest count equals the sum of partition records
- zero discarded findings
- zero database files
- no report limit or monolithic report timeout
- no SQLite in the single-pass BAM implementation
- preserved full BAM, BAI, BED, TSV, IGV session, and IGV batch outputs
