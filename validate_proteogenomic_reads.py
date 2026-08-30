#!/usr/bin/env python3
"""Read-level validation for PGTK proteogenomic findings.

Consumes the final variant/junction reports, Ensembl GTF, STAR coordinate-sorted
BAMs, and Arriba tables. Uses samtools for indexed extraction. Produces
consolidated event/read tables, compact IGV-ready BAMs, BED/BEDPE markers, and
an evidence report. No analyst-defined acceptance threshold is applied.
"""
import argparse, csv, gzip, json, re, subprocess, sys
import pysam
from collections import Counter, defaultdict
from pathlib import Path
from variant_read_evidence import classify_sam_fields

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
ATTR_RE = re.compile(r'(\S+) "([^"]*)";')
JUNCTION_RE = re.compile(r"([^:]+):(\d+)-(\d+)\(([+-])\)")


def run(cmd, capture=False):
    p = subprocess.run(cmd, text=True, capture_output=capture)
    if p.returncode:
        if capture and p.stderr: sys.stderr.write(p.stderr)
        raise SystemExit(f"command failed ({p.returncode}): {' '.join(map(str,cmd))}")
    return p.stdout if capture else ""


def open_text(path):
    p=Path(path)
    return gzip.open(p,"rt",encoding="utf-8",errors="replace") if p.suffix==".gz" else p.open(encoding="utf-8",errors="replace")


def read_tsv(path):
    with open_text(path) as h: return list(csv.DictReader(h,delimiter="\t"))


def write_tsv(path, rows, fields):
    with Path(path).open("w",encoding="utf-8",newline="") as h:
        w=csv.DictWriter(h,fieldnames=fields,delimiter="\t",lineterminator="\n",extrasaction="ignore"); w.writeheader(); w.writerows(rows)


def sample_from_path(path):
    sample = Path(path).name.split(".", 1)[0].strip()
    if not sample: raise SystemExit(f"cannot infer sample from {path}")
    return sample


def classify(ref,alt):
    if len(ref)==1 and len(alt)==1: return "SNV"
    if len(ref)==len(alt): return "MNV" if len(ref)>1 else "SNV"
    if len(ref)<len(alt): return "insertion" if alt.startswith(ref) else "complex_allele"
    if len(ref)>len(alt): return "deletion" if ref.startswith(alt) else "complex_allele"
    return "complex_allele"


def contigs(bam):
    return {x.split("\t",1)[0] for x in pysam.idxstats(str(bam)).splitlines() if x and not x.startswith("*")}


def match_contig(chrom, available):
    for x in (chrom, chrom[3:] if chrom.startswith("chr") else "chr"+chrom):
        if x in available:return x
    return None


def cigar_observation(seq,qual,cigar,start,pos,ref,alt):
    operations=[(int(n),op) for n,op in CIGAR_RE.findall(cigar)]
    # VCF indels are anchored at POS; CIGAR I/D normally begins after that anchor.
    rp=start; qp=0
    for n,op in operations:
        if op in "M=X": rp+=n; qp+=n
        elif op=="I":
            if len(alt)>len(ref) and rp-1==pos:
                return "+"+seq[qp:qp+n],None,qp+1,"insertion"
            qp+=n
        elif op=="D":
            if len(ref)>len(alt) and rp-1==pos:
                return f"-{n}",None,qp+1,"deletion"
            rp+=n
        elif op=="N": rp+=n
        elif op=="S": qp+=n
    rp=start; qp=0
    for n,op in operations:
        if op in "M=X":
            if rp<=pos<rp+n:
                i=qp+pos-rp; observed_length=max(1,len(ref)) if len(ref)==len(alt) else 1
                base=seq[i:i+observed_length]
                bq=(ord(qual[i])-33) if qual!="*" and i<len(qual) else None
                return base,bq,i+1,"base"
            rp+=n; qp+=n
        elif op in "DN":
            if rp<=pos<rp+n:return (f"-{n}" if op=="D" else f"N{n}"),None,qp+1,("deletion" if op=="D" else "splice_skip")
            rp+=n
        elif op in "IS": qp+=n
    return "",None,None,"none"

def support_class(obs,kind,ref,alt):
    if kind=="SNV":
        return "ALT" if obs.upper()==alt.upper() else "REF" if obs.upper()==ref.upper() else "OTHER"
    if kind=="MNV":
        return "ALT" if obs.upper().startswith(alt.upper()) else "REF" if obs.upper().startswith(ref.upper()) else "OTHER"
    if kind=="insertion":
        inserted=alt[len(ref):] if alt.startswith(ref) else alt
        return "ALT" if obs.upper()==("+"+inserted).upper() else "REF" if obs.upper().startswith(ref[0].upper()) else "OTHER"
    if kind=="deletion":
        deleted=len(ref)-len(alt)
        return "ALT" if obs==f"-{deleted}" else "REF" if obs.upper().startswith(ref[0].upper()) else "OTHER"
    return "OBSERVED"


def load_gtf_features(gtf, wanted):
    out=defaultdict(list)
    with open_text(gtf) as h:
        for line in h:
            if line.startswith("#"):continue
            f=line.rstrip().split("\t")
            if len(f)!=9:continue
            a=dict(ATTR_RE.findall(f[8])); gene=a.get("gene_name","")
            if gene in wanted:
                out[gene].append((f[0],int(f[3]),int(f[4]),f[2],f[6],a.get("transcript_id",""),a.get("exon_number","")))
    return out


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--variants",required=True); ap.add_argument("--junctions",required=True)
    ap.add_argument("--splice-detail",required=True); ap.add_argument("--arriba",nargs="*",default=[])
    ap.add_argument("--bam",action="append",required=True,metavar="SAMPLE=PATH")
    ap.add_argument("--gtf",required=True); ap.add_argument("--genome",required=True)
    ap.add_argument("--output-prefix",default="proteogenomic_read_validation")
    ap.add_argument("--padding",type=int,default=150); ap.add_argument("--min-mapping-quality",type=int,default=20); ap.add_argument("--min-base-quality",type=int,default=20)
    args=ap.parse_args()
    if not getattr(pysam, "__samtools_version__", ""):
        raise SystemExit("pysam does not expose an HTSlib/samtools runtime")
    prefix=Path(args.output_prefix); outdir=prefix.parent; outdir.mkdir(parents=True,exist_ok=True)
    variants=read_tsv(args.variants); junctions=read_tsv(args.junctions)
    selected=[r for r in variants if r.get("Altered-residue peptides")]
    genes={g for r in selected for g in r.get("Genes","").split(";") if g}
    gtf=load_gtf_features(args.gtf,genes)
    bam_map={x.split("=",1)[0]:Path(x.split("=",1)[1]) for x in args.bam}
    available={s:contigs(b) for s,b in bam_map.items()}

    event_rows=[]; read_rows=[]; bed_rows=[]
    for number,r in enumerate(selected,1):
        chrom=r["Chromosome"]; pos=int(r["Position"]); ref=r["REF"]; alt=r["ALT"]; kind=classify(ref,alt)
        eid=f"VAR{number:06d}_{r['Sample']}_{chrom}_{pos}_{ref}_{alt}".replace("/","_")
        for sample,bam in sorted(bam_map.items()):
            contig=match_contig(chrom,available[sample])
            if not contig: continue
            end=pos+max(len(ref),1)-1; region=f"{contig}:{max(1,pos-args.padding)}-{end+args.padding}"
            text=pysam.view(str(bam), region)
            counts=Counter(); mapqs=[]; bqs=[]; strands=Counter(); mates=Counter(); names=[]
            for line in text.splitlines():
                f=line.split("\t");
                if len(f)<11:continue
                flag=int(f[1]); start=int(f[3]); classification,otype,obs,qtext=classify_sam_fields(f,pos,ref,alt,args.min_mapping_quality,args.min_base_quality)
                support="ALT" if classification=="EXACT_ALT" else "REF" if classification=="CLEAN_REFERENCE" else "OTHER"; bq=int(qtext.split(";",1)[0]) if qtext and qtext.split(";",1)[0].isdigit() else None; rpos=""; counts[support]+=1; mapqs.append(int(f[4]));
                if bq is not None:bqs.append(bq)
                strand="-" if flag&16 else "+"; strands[(support,strand)]+=1
                mate="R1" if flag&64 else "R2" if flag&128 else "unpaired"; mates[(support,mate)]+=1
                names.append(f[0])
                read_rows.append({"Event ID":eid,"RNA event sample":r["Sample"],"BAM sample":sample,"Read name":f[0],"Contig":contig,"Alignment start":start,"MAPQ":f[4],"Flag":flag,"CIGAR":f[5],"Strand":strand,"Mate":mate,"Observation":obs,"Observation type":otype,"Support":support,"Base quality":"" if bq is None else bq,"Read position":"" if rpos is None else rpos})
            info=counts["REF"]+counts["ALT"]
            event_rows.append({"Event ID":eid,"RNA event sample":r["Sample"],"BAM sample":sample,"Genes":r.get("Genes",""),"Variant":r["Variant"],"Type":kind,"HGVSc":r.get("HGVSc",""),"HGVSp":r.get("HGVSp",""),"Peptides":r.get("Altered-residue peptides",""),"Contig":contig,"Start":pos,"End":end,"REF reads":counts["REF"],"ALT reads":counts["ALT"],"Other reads":sum(counts.values())-info,"ALT fraction":f"{counts['ALT']/info:.6f}" if info else "NA","ALT forward":strands[("ALT","+")],"ALT reverse":strands[("ALT","-")],"ALT R1":mates[("ALT","R1")],"ALT R2":mates[("ALT","R2")],"Mean MAPQ":f"{sum(mapqs)/len(mapqs):.2f}" if mapqs else "NA","Mean base quality":f"{sum(bqs)/len(bqs):.2f}" if bqs else "NA","Overlapping Ensembl111 features":sum(1 for g in r.get("Genes","").split(";") for x in gtf.get(g,[]) if x[1]<=pos<=x[2])})
        bed_rows.append((chrom,max(0,pos-1),pos+max(len(ref),1)-1,eid))

    event_fields=["Event ID","RNA event sample","BAM sample","Genes","Variant","Type","HGVSc","HGVSp","Peptides","Contig","Start","End","REF reads","ALT reads","Other reads","ALT fraction","ALT forward","ALT reverse","ALT R1","ALT R2","Mean MAPQ","Mean base quality","Overlapping Ensembl111 features"]
    read_fields=["Event ID","RNA event sample","BAM sample","Read name","Contig","Alignment start","MAPQ","Flag","CIGAR","Strand","Mate","Observation","Observation type","Support","Base quality","Read position"]
    write_tsv(str(prefix)+".events.tsv",event_rows,event_fields); write_tsv(str(prefix)+".reads.tsv",read_rows,read_fields)
    with open(str(prefix)+".variants.bed","w") as h:
        for x in bed_rows:h.write("\t".join(map(str,x))+"\n")

    # compact union BAMs for all altered-residue loci
    region_bed=Path(str(prefix)+".extraction_regions.bed")
    with region_bed.open("w") as h:
        seen=set()
        for chrom,start,end,eid in bed_rows:
            key=(chrom,max(0,start-args.padding),end+args.padding)
            if key not in seen: h.write(f"{key[0]}\t{key[1]}\t{key[2]}\n"); seen.add(key)
    bam_outputs=[]
    for sample,bam in sorted(bam_map.items()):
        ob=Path(f"{prefix}.{sample}.bam")
        pysam.view("-bh", "-L", str(region_bed), "-o", str(ob), str(bam), catch_stdout=False)
        pysam.index(str(ob)); bam_outputs.append(ob)

    # fusion BEDPE and consolidated Arriba table
    fusion_rows=[]
    with open(str(prefix)+".fusions.bedpe","w") as bedpe:
        for path in args.arriba:
            sample=sample_from_path(path)
            for r in read_tsv(path):
                b1=r.get("breakpoint1",""); b2=r.get("breakpoint2","")
                if ":" not in b1 or ":" not in b2:continue
                c1,p1=b1.rsplit(":",1); c2,p2=b2.rsplit(":",1)
                fid=f"FUSION_{sample}_{r.get('#gene1','')}_{r.get('gene2','')}_{c1}_{p1}_{c2}_{p2}"
                bedpe.write(f"{c1}\t{int(p1)-1}\t{p1}\t{c2}\t{int(p2)-1}\t{p2}\t{fid}\n")
                fusion_rows.append({"Event ID":fid,"Sample":sample,"Gene 1":r.get("#gene1",""),"Gene 2":r.get("gene2",""),"Breakpoint 1":b1,"Breakpoint 2":b2,"Type":r.get("type",""),"Confidence":r.get("confidence",""),"Split reads 1":r.get("split_reads1",""),"Split reads 2":r.get("split_reads2",""),"Discordant mates":r.get("discordant_mates",""),"Reading frame":r.get("reading_frame",""),"Filters":r.get("filters","")})
    ff=["Event ID","Sample","Gene 1","Gene 2","Breakpoint 1","Breakpoint 2","Type","Confidence","Split reads 1","Split reads 2","Discordant mates","Reading frame","Filters"]
    write_tsv(str(prefix)+".fusions.tsv",fusion_rows,ff)

    # junction table and BED markers
    with open(str(prefix)+".junctions.bed","w") as h:
        for i,r in enumerate(junctions,1):
            m=JUNCTION_RE.match(r.get("Genomic junction",""))
            if m:h.write(f"{m.group(1)}\t{int(m.group(2))-1}\t{m.group(3)}\tJUNCTION{i:04d}_{r.get('Sequence','')}\n")

    igv=Path(str(prefix)+".igv.batch.txt")
    with igv.open("w") as h:
        h.write("new\ngenome hg38\n")
        for ob in bam_outputs:h.write(f"load {ob.resolve()}\n")
        h.write(f"load {Path(str(prefix)+'.variants.bed').resolve()}\n")
        h.write(f"load {Path(str(prefix)+'.fusions.bedpe').resolve()}\n")
        h.write(f"load {Path(str(prefix)+'.junctions.bed').resolve()}\n")

    summary=Counter(r["Type"] for r in event_rows if r["BAM sample"]==r["RNA event sample"])
    with open(str(prefix)+".summary.txt","w") as h:
        h.write(f"Altered-residue genomic events: {len(selected)}\n")
        h.write(f"Read-level sample-event comparisons: {len(event_rows)}\n")
        h.write(f"Read observations: {len(read_rows)}\n")
        h.write(f"Accepted Arriba events: {len(fusion_rows)}\n")
        h.write(f"Translated junction findings: {len(junctions)}\n")
        for k in sorted(summary):h.write(f"{k} events: {summary[k]}\n")
    with open(str(prefix)+".report.md","w") as h:
        h.write("# Proteogenomic read validation\n\n")
        h.write("Read-level evidence is reported without analyst-defined acceptance thresholds. RNA observations do not establish somatic or genomic structural-variant status.\n\n")
        h.write(Path(str(prefix)+".summary.txt").read_text().replace("\n","  \n"))
        h.write("\n## Outputs\n\nEvents, reads, compact BAMs, BED/BEDPE markers, fusion summaries, junction markers, and an IGV batch file are included.\n")

if __name__=="__main__":main()
