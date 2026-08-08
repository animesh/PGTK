#!/usr/bin/env python3
import re
from pathlib import Path
root=Path(__file__).resolve().parent
main=(root/'main.nf').read_text()
config=(root/'nextflow.config').read_text()
slurm=(root/'scratch.slurm').read_text()
processes=re.findall(r'^process\s+(\w+)\s*\{',main,re.M)
robust=config.split('    robust {',1)[1].rsplit('\n    }\n}',1)[0]
selectors=re.findall(r'withName:\s*(\w+)\s*\{',robust)
assert len(processes)==53 and len(set(processes))==53
assert set(processes)==set(selectors),(set(processes)-set(selectors),set(selectors)-set(processes))
assert "queue = {" in config
assert 'task.cpus > pgtkEffectiveNormalCpuThreshold || task.memory > pgtkEffectiveNormalMemoryThresholdGb.GB' in config
assert "requiredEnv('PGTK_NORMAL_PARTITION')" in config
assert "requiredEnv('PGTK_BIGMEM_PARTITION')" in config
assert "requiredEnv('PGTK_MAX_CPUS')" in config
assert "requiredEnv('PGTK_MAX_MEMORY_GB')" in config
assert 'Math.min(pgtkMaxCpus' in config
assert 'value > pgtkEffectiveMaxMemoryGb.GB ? pgtkEffectiveMaxMemoryGb.GB : value' in config
assert 'pgtkAbsoluteNormalCpuThreshold = 20' in config
assert 'pgtkAbsoluteNormalMemoryThresholdGb = 160' in config
assert 'pgtkAbsoluteMaxCpus = 32' in config
assert 'pgtkAbsoluteMaxMemoryGb = 512' in config
assert 'maxRetries = 2' in config
assert "task.exitStatus in [137, 140, 143] ? 'retry' : 'terminate'" in config
assert not re.search(r"queue\s*=\s*'(normal|bigmem)'",config)
assert not re.search(r"queue\s+'(normal|bigmem)'",main)
assert '#SBATCH --account=nn9036k' in slurm
assert '#SBATCH --partition=normal' in slurm
for token in ['--project-dir','--nextflow','--python','--apptainer','--work-dir','--tmp-root','--account','--normal-partition','--bigmem-partition','--normal-cpu-threshold','--normal-memory-threshold-gb','--max-cpus','--max-memory-gb','--queue-size','--submit-rate-limit','--slurm-log-template']:
    assert token in slurm,token
for required in [
    'PROJECT_DIR=/cluster/projects/nn9036k/scrbkup/pgtk',
    'NEXTFLOW=/cluster/home/ash022/bin/nextflow',
    'ACCOUNT=nn9036k',
    'NORMAL_PARTITION=normal',
    'BIGMEM_PARTITION=bigmem',
    'JAVA_MODULE=Java/21',
]:
    assert required in slurm, required
assert 'PIPELINE_ARGS=("$@")' in slurm
assert 'trap finalize EXIT' in slurm
assert '--host_python "$HOST_PYTHON"' in slurm
assert '--host-python "$HOST_PYTHON"' in slurm
assert '--apptainer "$APPTAINER"' in slurm
assert '"${params.host_python}" ${bundle_script}' in main

def allocation(base_cpu,base_mem,attempt,cpu_scales):
    multiplier=1 << (attempt-1)
    cpu=min(32,base_cpu*multiplier) if cpu_scales else base_cpu
    mem=min(512,base_mem*multiplier)
    queue='bigmem' if cpu>20 or mem>160 else 'normal'
    return cpu,mem,queue
assert allocation(2,48,1,False)==(2,48,'normal')
assert allocation(2,48,2,False)==(2,96,'normal')
assert allocation(2,48,3,False)==(2,192,'bigmem')
assert allocation(32,64,1,True)==(32,64,'bigmem')
assert allocation(32,64,3,True)==(32,256,'bigmem')
assert allocation(3,6,3,True)==(12,24,'normal')
assert allocation(1,128,2,False)==(1,256,'bigmem')
assert main.count('def javaHeapGb = Math.max(1, Math.floor(task.memory.toGiga() * 0.80) as int)') == 6
assert '-Xmx40g' not in main and '-Xmx18g' not in main and '-Xmx6g' not in main and '-Xmx4g' not in main
assert main.count('-Xmx${javaHeapGb}g') == 9
assert 'ParallelGCThreads=${javaGcThreads}' in main
for process in ['MARK_DUPLICATES','SPLIT_N_CIGAR','PREPARE_HAPLOTYPE_INTERVALS','GATHER_HAPLOTYPE_GVCF','GENOTYPE_FILTER']:
    block=config.split(f'withName: {process} {{',1)[1].split('\n            }',1)[0]
    assert 'cpus = { Math.min' not in block, process
for process in ['HAPLOTYPE_CALLER','STAR_ALIGN','STAR_INDEX','SORT_INDEX_BAM']:
    block=config.split(f'withName: {process} {{',1)[1].split('\n            }',1)[0]
    assert 'cpus = { Math.min' in block, process
print('PASS: dynamic partition routing, retry scaling, caps and CLI-only wrapper configuration')
python_script_variables = ('${validator}', '${validation_script}', '${bundle_script}')
for match in re.finditer(r'^process\s+(\w+)\s*\{', main, re.M):
    process = match.group(1)
    following = re.search(r'^process\s+\w+\s*\{', main[match.end():], re.M)
    end = match.end() + following.start() if following else len(main)
    block = main[match.start():end]
    if 'samtools:1.21' in block and any(variable in block for variable in python_script_variables):
        assert '${params.host_python}' in block, f'{process} uses a Python script in the samtools container without host Python'
print('PASS: samtools-container Python consumers use configured host Python')
