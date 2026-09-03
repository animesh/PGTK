#!/usr/bin/env python3
import ast
import re
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def imported_call_contracts():
    modules={p.stem:p for p in ROOT.glob('*.py')}
    defs={}
    for name,path in modules.items():
        tree=ast.parse(path.read_text(),filename=str(path))
        defs[name]={n.name:n for n in tree.body if isinstance(n,ast.FunctionDef)}
    errors=[]
    for path in modules.values():
        tree=ast.parse(path.read_text(),filename=str(path)); imports={}
        for n in ast.walk(tree):
            if isinstance(n,ast.ImportFrom) and n.module in modules:
                for item in n.names: imports[item.asname or item.name]=(n.module,item.name)
        for n in ast.walk(tree):
            if isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id in imports:
                mod,fun=imports[n.func.id]; d=defs.get(mod,{}).get(fun)
                if d is None: errors.append(f'{path.name}:{n.lineno}: missing {mod}.{fun}'); continue
                positional=d.args.posonlyargs+d.args.args
                minimum=len(positional)-len(d.args.defaults)
                maximum=None if d.args.vararg else len(positional)
                supplied=len(n.args)
                if supplied<minimum or (maximum is not None and supplied>maximum):
                    errors.append(f'{path.name}:{n.lineno}: {mod}.{fun} accepts {minimum}..{maximum}; got {supplied}')
    assert not errors,'\n'.join(errors)

def production_classifier_contract():
    import sys
    sys.path.insert(0,str(ROOT))
    from build_finding_igv_reviews import classify
    fields=['r','0','1','100','60','10M','*','0','0','ACGTACGTAC','IIIIIIIIII']
    result=classify(fields,101,'C','T',20,20)
    assert result[0]=='CLEAN_REFERENCE',result

def nextflow_python_cli_contracts():
    source=(ROOT/'main.nf').read_text()
    options={}
    for path in ROOT.glob('*.py'):
        tree=ast.parse(path.read_text(),filename=str(path)); found=set()
        for n in ast.walk(tree):
            if isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and n.func.attr=='add_argument':
                for a in n.args:
                    if isinstance(a,ast.Constant) and isinstance(a.value,str) and a.value.startswith('--'): found.add(a.value)
        options[path.name]=found
    mappings=re.findall(r'(\w+)\s*=\s*file\("\$\{projectDir\}/([^"}]+\.py)"',source)
    variables=dict(mappings)
    errors=[]
    for match in re.finditer(r'python3?\s+\$\{(\w+)\}',source):
        variable=match.group(1); script=variables.get(variable)
        if not script: continue
        end=source.find('"""',match.end())
        command=source[match.start():end if end!=-1 else match.end()+2000]
        command=re.sub(r'\$\{[^}]+\}','',command)
        supplied=set(re.findall(r'--[A-Za-z][A-Za-z0-9_-]*',command))
        unknown=supplied-options.get(script,set())
        if unknown: errors.append(f'{script}: unsupported options {sorted(unknown)}')
    # Dynamic argparse declarations are covered by executable --help tests in runtime preflight.
    # This static check is intentionally limited to scripts with literal declarations.
    literal_errors=[e for e in errors if not any(name in e for name in ('compare_progression_pair.py','merge_progression_biology.py','build_compact_multiqc_content.py'))]
    assert not literal_errors,'\n'.join(literal_errors)

imported_call_contracts()
production_classifier_contract()


def samplesheet_validator_contract():
    import csv
    import subprocess
    import sys
    import tempfile
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        source = root / "samples.csv"
        validated = root / "validated.csv"
        report = root / "report.tsv"
        source.write_text("sample,srr,TK,Group,baseline\nA,SRR1,TK1,G,true\nB,SRR2,TK1,G,false\n", encoding="utf-8")
        subprocess.run([sys.executable, str(ROOT / "validate_samplesheet_design.py"), "--input", str(source), "--validated-samplesheet", str(validated), "--report", str(report)], check=True)
        with report.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        assert rows[0]["SubtractionStatus"] == "ENABLED"
        source.write_text("sample,srr,TK,Group,baseline\nA,SRR1,TK1,G,true\nB,SRR2,TK1,G,true\n", encoding="utf-8")
        failed = subprocess.run([sys.executable, str(ROOT / "validate_samplesheet_design.py"), "--input", str(source), "--validated-samplesheet", str(validated), "--report", str(report)], capture_output=True, text=True)
        assert failed.returncode != 0
        assert "multiple baselines" in failed.stderr



def complete_validator_process_name_contract():
    import importlib.util
    specification = importlib.util.spec_from_file_location("complete_validator", ROOT / "validate_pgtk_results_complete.py")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    cases = {
        "HAPLOTYPE_CALLER (TK12:0000-scattered)": "HAPLOTYPE_CALLER",
        "PGTK:HAPLOTYPE_CALLER (TK12:0000-scattered)": "HAPLOTYPE_CALLER",
        "ANALYZE_PROGRESSION_SAMPLE (TK14:GO_progression)": "ANALYZE_PROGRESSION_SAMPLE",
        "SRA_TO_FASTQ (TK12:SRR31089074)": "SRA_TO_FASTQ",
        "STAR_INDEX (GRCh38_Ensembl111)": "STAR_INDEX",
        "VALIDATE_PUBLISHED_RESULTS (postrun_validation:19710852)": "VALIDATE_PUBLISHED_RESULTS",
    }
    for trace_name, expected in cases.items():
        observed = module.process_name(trace_name)
        assert observed == expected, (trace_name, observed, expected)


def complete_validator_pipeline_mode_contract():
    import json
    import subprocess
    import sys
    import tempfile
    validator = ROOT / "validate_pgtk_results_complete.py"
    with tempfile.TemporaryDirectory() as directory:
        temporary = Path(directory)
        results = temporary / "results"
        output_root = temporary / "output"
        results.mkdir(); output_root.mkdir()
        job_id = "fixture"
        (results / f"pipeline_trace-{job_id}.tsv").write_text("name\tstatus\nHAPLOTYPE_CALLER (A:0000-scattered)\tCACHED\nVALIDATE_PUBLISHED_RESULTS (postrun_validation:fixture)\tRUNNING\n", encoding="utf-8")
        deep = output_root / f"PGTK-deep-audit-{job_id}"
        deep.mkdir()
        (deep / "summary.json").write_text(json.dumps({"audit_workers": 32, "display_events_selected": 0, "findings": 0}), encoding="utf-8")
        (deep / "issues.tsv").write_text("level\tcode\tmessage\tpath\tdetails\n", encoding="utf-8")
        completed = subprocess.run([sys.executable, str(validator), "--pipeline-mode", "--project-dir", str(ROOT), "--results-dir", str(results), "--job-id", job_id, "--workers", "32", "--max-events", "0", "--output-root", str(output_root), "--pysam-image", str(validator), "--host-python", sys.executable, "--nextflow", "/bin/false", "--apptainer", "/bin/false", "--reuse-deep-audit"], cwd=ROOT, capture_output=True, text=True)
        assert completed.returncode == 2, completed.stdout + completed.stderr
        assert "Traceback" not in completed.stdout + completed.stderr
        assert "unknown=['0000-scattered)']" not in completed.stdout + completed.stderr
        report_dir = output_root / f"PGTK-complete-validation-{job_id}"
        for name in ("REPORT.md", "checks.tsv", "summary.json"):
            assert (report_dir / name).is_file(), name
        assert (output_root / f"PGTK-complete-validation-{job_id}.tar.gz").is_file()
        assert (output_root / f"PGTK-complete-validation-{job_id}.tar.gz.sha256").is_file()

samplesheet_validator_contract()
complete_validator_process_name_contract()
complete_validator_pipeline_mode_contract()
print('program interface contract tests: PASS')
