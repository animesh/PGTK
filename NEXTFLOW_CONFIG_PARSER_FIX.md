# Nextflow configuration parser fix

Nextflow 25.04.4 rejected top-level `def` variable declarations mixed with configuration statements.

The fix removes only the 18 top-level `def` keywords. It preserves runtime environment resolution, partition routing, CPU and memory limits, retries, queue settings, all 64 process selectors, ranked GO safeguards, and GO FDR default 0.1.

Closure-local variables remain unchanged.
