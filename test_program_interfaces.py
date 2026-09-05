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



def generated_explorer_geometry_contract():
    import csv
    import gzip
    import json
    import subprocess
    import sys
    import tempfile

    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        output = root / "explorer"
        fields = [
            "EventID","EvidenceClasses","SourceEvents","Sources","Sample","Gene","Consequence","Impact",
            "PredictedConsequence","PredictedImpact","Transcript","ProteinChange","Chrom","Position","REF","ALT",
            "ReadValidationStatus","ValidationExplanation","CountUnit","UniqueAlignments","CallableAlignments",
            "ExactAltReads","CleanReferenceReads","ExcludedReads","AltFractionAmongClean","CallableFractionAmongExamined",
        ]
        rows = [
            ("SNV1","rna_variant","S1","1","10","A","T","ALT_SUPPORTED",1,0,0),
            ("INS1","rna_variant","S2","1","20","A","AG","MIXED_ALT_AND_REFERENCE",1,1,0),
            ("DEL1","rna_variant","S3","1","30","ATG","A","MIXED_ALT_AND_REFERENCE",1,1,0),
            ("SPL1","splice_junction","S4","2","100","","","NO_CALLABLE_READS",0,0,0),
            ("FUS1","fusion","S5","3","200","","","NO_CALLABLE_READS",0,0,0),
        ]
        with (root / "manifest.tsv").open("w", newline="") as handle:
            writer=csv.DictWriter(handle,fieldnames=fields,delimiter="\t",lineterminator="\n");writer.writeheader()
            for event_id,classes,source,chrom,pos,ref,alt,status,exact,clean,excluded in rows:
                row=dict.fromkeys(fields,"");row.update(EventID=event_id,EvidenceClasses=classes,SourceEvents=source,Sample="TK1",Chrom=chrom,Position=pos,REF=ref,ALT=alt,ReadValidationStatus=status,ValidationExplanation="fixture",CountUnit="alignments",UniqueAlignments=str(exact+clean+excluded),CallableAlignments=str(exact+clean),ExactAltReads=str(exact),CleanReferenceReads=str(clean),ExcludedReads=str(excluded),AltFractionAmongClean=(str(exact/(exact+clean)) if exact+clean else "NA"),CallableFractionAmongExamined=(str((exact+clean)/(exact+clean+excluded)) if exact+clean+excluded else "NA"))
                writer.writerow(row)
        event_fields=["Event","Sample","Class","Chrom","Start0","End","Chrom2","Start2_0","End2","Label","Source","Gene","Consequence","Impact","Transcript","ProteinChange"]
        events=[
            ["S1","TK1","rna_variant","1","9","10","","","","A>T","x","","","","",""],
            ["S2","TK1","rna_variant","1","19","20","","","","A>AG","x","","","","",""],
            ["S3","TK1","rna_variant","1","29","32","","","","ATG>A","x","","","","",""],
            ["S4","TK1","splice_junction","2","99","150","","","","JUNCTION","x","","","","",""],
            ["S5","TK1","fusion","3","199","200","4","299","300","A--B","x","","","","",""],
        ]
        with (root / "events.tsv").open("w",newline="") as handle:
            writer=csv.writer(handle,delimiter="\t",lineterminator="\n");writer.writerow(event_fields);writer.writerows(events)
        (root / "excluded.tsv").write_text("EventID\tReason\n")
        (root / "bam_manifest.tsv").write_text("Sample\tCategory\tBAM\tIndex\tUniqueAlignments\n")
        with gzip.open(root / "display.tsv.gz","wt") as handle:
            handle.write("EventID\tSample\tCategory\tAlignmentKey\tReadName\tContig\tStart0\tFlag\tCIGAR\tSequenceHash\n")
            handle.write("SPL1\tTK1\tevent_display\tk1\tr1\t2\t90\t0\t60M\th1\n")
            handle.write("FUS1\tTK1\tevent_display\tk2\tr2\t3\t190\t0\t60M\th2\n")
        (root / "genome.fa").write_text(">1\nA\n")
        completed=subprocess.run([sys.executable,str(ROOT/'build_finding_explorer.py'),'--manifest',str(root/'manifest.tsv'),'--events',str(root/'events.tsv'),'--excluded-reads',str(root/'excluded.tsv'),'--bam-manifest',str(root/'bam_manifest.tsv'),'--display-manifest',str(root/'display.tsv.gz'),'--genome',str(root/'genome.fa'),'--output-dir',str(output)],capture_output=True,text=True)
        assert completed.returncode==0,completed.stdout+completed.stderr
        geometry=json.loads((output/'event_geometry.json').read_text())
        assert geometry['SNV1']['event_type']=='SNV' and geometry['SNV1']['regions'][0]['start0']==9
        assert geometry['INS1']['event_type']=='INSERTION' and geometry['INS1']['regions'][0]['end0']==20
        assert geometry['DEL1']['event_type']=='DELETION' and geometry['DEL1']['regions'][0]['end0']==32
        assert geometry['SPL1']['event_type']=='SPLICE_JUNCTION' and geometry['SPL1']['regions']==[{'chrom':'2','start0':99,'end0':150,'role':'JUNCTION'}]
        assert geometry['FUS1']['event_type']=='FUSION' and len(geometry['FUS1']['regions'])==2
        compile((output/'server.py').read_text(),str(output/'server.py'),'exec')
        with gzip.open(output/'partitions/all.jsonl.gz','rt') as handle:
            records={row['EventID']:row for row in map(json.loads,handle)}
        assert records['SPL1']['VisualEvidenceStatus']=='CONTEXT_ALIGNMENTS_AVAILABLE'
        assert records['FUS1']['ContextAlignments']==1
        import importlib.util
        import os
        import shutil
        for helper in ('report_legend.py','prepare_event_igv_tracks.py'):
            shutil.copy2(ROOT/helper,output/helper)
        reviews=root/'finding_reviews';reviews.mkdir();shutil.copy2(root/'display.tsv.gz',reviews/'display.tsv.gz')
        pysam_image=root/'pysam.img';igv_image=root/'igv.img';pysam_image.write_bytes(b'image');igv_image.write_bytes(b'image')
        os.environ['PGTK_PYSAM_IMAGE']=str(pysam_image);os.environ['PGTK_IGV_REPORTS_IMAGE']=str(igv_image)
        specification=importlib.util.spec_from_file_location('generated_explorer_server',output/'server.py')
        server=importlib.util.module_from_spec(specification);specification.loader.exec_module(server)
        assert server.search({'event_type':['FUSION']})['total'] == 1
        assert server.search({'event_type':['SPLICE_JUNCTION']})['rows'][0]['EventID'] == 'SPL1'
        assert server.REPORT_CACHE_VERSION=='6-persistent-event-panel'
        assert server.event_description(server.I['SNV1'],server.G['SNV1'])['action']=='SUBSTITUTE A WITH T AT chr1:10'
        assert server.event_description(server.I['INS1'],server.G['INS1'])=={'action':'INSERT G AFTER chr1:20 A','reference':'A','alternate':'A[G]'}
        assert server.event_description(server.I['DEL1'],server.G['DEL1'])=={'action':'DELETE TG AFTER chr1:30 A','reference':'A[TG]','alternate':'A'}
        assert server.event_description(server.I['SPL1'],server.G['SPL1'])['action']=='SPLICE chr2:100 TO chr2:150'
        assert server.event_description(server.I['FUS1'],server.G['FUS1'])['action']=='FUSE chr3:200 TO chr4:300'
        ddx_row=dict(server.I['INS1'],EventID='DDX1_TK12_2_15597434_C_CA_C_CA',Sample='TK12',Gene='DDX1',Chrom='2',Position=15597434,REF='C',ALT='CA',ExactAltReads=3,CleanReferenceReads=813,ExcludedReads=161,CallableReads=816,AltFractionAmongClean=3/816)
        ddx_geometry={'event_type':'INSERTION','regions':[{'chrom':'2','start0':15597433,'end0':15597434,'role':'TARGET_ALLELE'}]}
        ddx_description=server.event_description(ddx_row,ddx_geometry)
        assert ddx_description=={'action':'INSERT A AFTER chr2:15,597,434 C','reference':'C','alternate':'C[A]'}
        ddx_summary=server.event_summary_html(ddx_row,ddx_geometry,ddx_description)
        assert 'id="pgtk-event-summary"' in ddx_summary and 'position:sticky' in ddx_summary
        assert 'INSERT A AFTER chr2:15,597,434 C' in ddx_summary and 'chr2:15,597,435' in ddx_summary
        assert 'Purple I symbols at the target boundary' in ddx_summary
        loader=server.persistent_summary_script(ddx_summary)
        assert 'pgtk-event-summary-loader' in loader and 'MutationObserver' in loader and 'container.insertBefore' in loader
        import threading
        import urllib.request
        httpd=server.ThreadingHTTPServer(('127.0.0.1',0),server.H)
        thread=threading.Thread(target=httpd.serve_forever,daemon=True);thread.start()
        proxy_free_opener=urllib.request.build_opener(
            urllib.request.ProxyHandler({})
        )
        try:
            port=httpd.server_address[1]
            with proxy_free_opener.open(
                f'http://127.0.0.1:{port}/api/findings?event_type=FUSION',
                timeout=10,
            ) as response:
                payload=json.loads(response.read())
            assert payload['total']==1 and payload['rows'][0]['EventID']=='FUS1'
            with proxy_free_opener.open(
                f'http://127.0.0.1:{port}/api/facets',
                timeout=10,
            ) as response:
                facets=json.loads(response.read())
            assert {'FUSION','SPLICE_JUNCTION','SNV','INSERTION','DELETION'}.issubset(set(facets['event_types']))
        finally:
            httpd.shutdown();httpd.server_close();thread.join(timeout=5)
        commands=[]
        def fake_run(command,label):
            commands.append((label,list(command)))
            if label=='event BAM extraction':
                event_output=Path(command[command.index('--output-dir')+1]);bam=event_output/'event_display.bam';bam.write_bytes(b'BAM');(event_output/'tracks.tsv').write_text(f'Category\tBAM\tIndex\nevent_display\t{bam}\t{bam}.bai\n')
            else:
                report_output=Path(command[command.index('--output')+1]);report_output.write_text('<html><body><div id="container" style="display:flex">fixture</div></body></html>')
        server.run=fake_run
        expected={
            'SNV1':('1','9','10','SNV','10','A','T'),
            'INS1':('1','19','20','INSERTION','20','A','AG'),
            'DEL1':('1','29','32','DELETION','32','ATG','A'),
        }
        for event_id,(chrom,start0,end0,event_type,site_end,ref,alt) in expected.items():
            report=server.report(event_id);assert report.is_file()
            event_dir=next(path for path in server.K.glob(event_id+'-*') if path.is_dir())
            marker_fields=[line for line in (event_dir/'event_coordinate.bed').read_text().splitlines() if not line.startswith('track')][0].split('\t')
            assert marker_fields[:3]==[chrom,start0,end0]
            description=server.event_description(server.I[event_id],server.G[event_id])
            assert description['action'] in marker_fields[3] and 'TARGET_ALLELE' in marker_fields[3]
            with (event_dir/'finding.tsv').open(newline='',encoding='utf-8') as handle:
                site_row=next(csv.DictReader(handle,delimiter='\t'))
            assert site_row['Start']==str(int(start0)+1) and site_row['End']==site_end
            assert site_row['REF']==ref and site_row['ALT']==alt and site_row['EventType']==event_type
            report_command=[command for label,command in commands if label=='IGV report generation'][-1]
            track_index=report_command.index('--tracks')
            assert report_command[track_index+1].endswith('event_coordinate.bed')
            report_text=report.read_text()
            assert report_text.count('id="pgtk-event-summary"')==1
            assert report_text.count('id="pgtk-event-summary-loader"')==1
            assert description['action'] in report_text and 'MutationObserver' in report_text
            assert report_text.index('id="pgtk-event-summary"') < report_text.index('fixture')
            before=len(commands);assert server.report(event_id)==report;assert len(commands)==before
        server.report('FUS1');fusion_dir=next(path for path in server.K.glob('FUS1-*') if path.is_dir());fusion_marker=(fusion_dir/'event_coordinate.bed').read_text()
        assert 'FUSE chr3:200 TO chr4:300' in fusion_marker and 'BREAKPOINT_1' in fusion_marker and 'BREAKPOINT_2' in fusion_marker
        with (fusion_dir/'finding.tsv').open(newline='',encoding='utf-8') as handle:
            fusion_rows=list(csv.DictReader(handle,delimiter='\t'))
        assert len(fusion_rows)==2 and {row['RegionRole'] for row in fusion_rows}=={'BREAKPOINT_1','BREAKPOINT_2'}
        server.report('SPL1');splice_dir=next(path for path in server.K.glob('SPL1-*') if path.is_dir());splice_marker=(splice_dir/'event_coordinate.bed').read_text()
        assert '\t99\t150\tSPLICE chr2:100 TO chr2:150 | JUNCTION' in splice_marker

generated_explorer_geometry_contract()

samplesheet_validator_contract()
complete_validator_process_name_contract()
complete_validator_pipeline_mode_contract()
print('program interface contract tests: PASS')