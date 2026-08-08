#!/usr/bin/env python3
import argparse, csv, gzip, hashlib
from collections import defaultdict
from pathlib import Path


def op(path):
    return gzip.open(path,'rt',encoding='utf-8',errors='replace') if str(path).endswith('.gz') else open(path,encoding='utf-8',errors='replace')


def parse_obo(path):
    terms={}; current=None; data_version='' 
    def save(term):
        if term and term.get('id') and not term.get('obsolete'):
            terms[term['id']]=term
    with op(path) as h:
        for raw in h:
            line=raw.rstrip('\n')
            if line.startswith('data-version:'): data_version=line.split(':',1)[1].strip()
            if line=='[Term]': save(current); current={'parents':set(),'alt_ids':set()}; continue
            if line.startswith('['): save(current); current=None; continue
            if current is None or ': ' not in line: continue
            k,v=line.split(': ',1)
            if k=='id': current['id']=v
            elif k=='name': current['name']=v
            elif k=='namespace': current['namespace']=v
            elif k=='alt_id': current['alt_ids'].add(v)
            elif k=='is_a': current['parents'].add(v.split()[0])
            elif k=='relationship' and v.startswith('part_of '): current['parents'].add(v.split()[1])
            elif k=='is_obsolete' and v=='true': current['obsolete']=True
    save(current)
    aliases={alt:go for go,t in terms.items() for alt in t['alt_ids']}
    return terms,aliases,data_version


def main():
    p=argparse.ArgumentParser(); p.add_argument('--obo',required=True); p.add_argument('--gaf',required=True); p.add_argument('--output-prefix',required=True); a=p.parse_args()
    terms,aliases,ontology_version=parse_obo(a.obo)
    direct=defaultdict(set); gaf_version=''; generated=''
    with op(a.gaf) as h:
        for line in h:
            if line.startswith('!gaf-version:'): gaf_version=line.split(':',1)[1].strip()
            if line.startswith('!generated-on:'): generated=line.split(':',1)[1].strip()
            if line.startswith('!'): continue
            f=line.rstrip('\n').split('\t')
            if len(f)<9 or 'NOT' in f[3].split('|'): continue
            gene=f[2].strip(); go=aliases.get(f[4],f[4]); aspect=f[8]
            if gene and go in terms and aspect in {'P','F','C'}: direct[gene].add(go)
    memo={}
    def ancestors(go,trail=frozenset()):
        if go in memo:return memo[go]
        if go in trail:return {go}
        out={go}
        for parent in terms.get(go,{}).get('parents',()):
            if parent in terms: out |= ancestors(parent,trail|{go})
        memo[go]=out; return out
    rows=[]
    for gene,gos in sorted(direct.items()):
        propagated=set().union(*(ancestors(go) for go in gos))
        for go in sorted(propagated):
            t=terms[go]; rows.append({'Gene':gene,'GO_ID':go,'GO_Name':t['name'],'Namespace':t['namespace']})
    out=a.output_prefix+'.mapping.tsv'
    with open(out,'w',encoding='utf-8',newline='') as h:
        w=csv.DictWriter(h,fieldnames=['Gene','GO_ID','GO_Name','Namespace'],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    defhash=lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
    Path(a.output_prefix+'.metadata.tsv').write_text('Field\tValue\nOntologyFile\t'+Path(a.obo).name+'\nOntologyVersion\t'+ontology_version+'\nOntologySHA256\t'+defhash(a.obo)+'\nAnnotationFile\t'+Path(a.gaf).name+'\nAnnotationSHA256\t'+defhash(a.gaf)+'\nGAFVersion\t'+gaf_version+'\nGeneratedOn\t'+generated+'\nOntologyTerms\t'+str(len(terms))+'\nAnnotatedGenes\t'+str(len(direct))+'\nPropagatedGeneTermPairs\t'+str(len(rows))+'\nPropagationRelations\tis_a,part_of\n',encoding='utf-8')
if __name__=='__main__':main()
