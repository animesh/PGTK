# Portable post-processing fix

Reviewed failure: production single-pass IGV computation completed for 157,482 findings, then the process exited because the Pysam container provides a `find` implementation without GNU `-printf`.

Changes:

- Replaced `find ... -printf` entry counting with POSIX-compatible `find ... | wc -l`.
- Replaced another diagnostic `find ... -printf` call in `GFFCOMPARE_NOVEL` with a portable shell loop and `wc -c`.
- Replaced GNU-specific `stat -c%s` in `GENERATE_PRIORITY_IGV_REPORTS` with portable `wc -c` redirection.
- Added regression assertions preventing reintroduction of `find -printf` and `stat -c` in `main.nf`.
- Preserved the validated single-pass BAM implementation, all 72 processes, resource configuration, samplesheet semantics, reporting branches, and test fixtures.
