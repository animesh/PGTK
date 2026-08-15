# Nextflow 26.04 strict configuration fix

The exact Saga diagnostic used Nextflow 26.04.6, where the strict parser is enabled by default.

The strict configuration language accepts configuration assignments, blocks and includes. It does not accept a top-level helper closure such as `requiredEnv`, and configuration settings cannot be referenced as ordinary variables.

This replacement uses the official `env('NAME')` configuration function directly inside configuration expressions. The wrapper remains responsible for validating that every required environment value is supplied before Nextflow starts.

Preserved behavior:

- dynamic normal versus bigmem routing
- runtime account and partition settings
- CPU cap of 32 and memory cap of 512 GB
- normal thresholds capped at 20 CPUs and 160 GB
- queue size and submission-rate settings
- three total attempts with resource-only retries
- all 64 process selectors
- ranked GO safeguards and FDR default 0.1
