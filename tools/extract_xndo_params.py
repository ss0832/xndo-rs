#!/usr/bin/env python3
"""Extract MNDO, MNDO/d, and ZINDO/S parameters from OpenMOPAC v23.2.5.

Usage:
  python tools/extract_xndo_params.py \
    --mndo /path/parameters_for_mndo_C.F90 \
    --mndod /path/parameters_for_mndod_C.F90 \
    --mndo-pairs /path/alpb_and_xfac_mndo.F90 \
    --mndod-pairs /path/alpb_and_xfac_mndod.F90 \
    --indo /path/parameters_for_INDO_C.F90
"""
from __future__ import annotations
import argparse, csv, re
from pathlib import Path

TAG = "v23.2.5"
REPO = "https://github.com/openmopac/mopac"
NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[DdEe][-+]?\d+)?"
SYMBOLS = [
    "H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca",
    "Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr",
    "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd",
    "Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg",
    "Tl","Pb","Bi","Po","At","Rn"
]

def norm(x:str)->str:
    return x.replace('D','e').replace('d','e')

def parse_data(path:Path):
    pat=re.compile(
        rf"^\s*data\s+([A-Za-z_]\w*)\s*\(\s*(\d+)\s*(?:,\s*(\d+)\s*)?\)\s*/\s*({NUM})\s*/",
        re.I,
    )
    out={}
    for line in path.read_text().splitlines():
        m=pat.match(line)
        if m:
            out[(m.group(1).lower(),int(m.group(2)),int(m.group(3)) if m.group(3) else None)]=norm(m.group(4))
    return out

def parse_pairs(path:Path):
    pat=re.compile(rf"^\s*(alpb|xfac)\s*\(\s*(\d+)\s*,\s*(\d+)\s*\)\s*=\s*({NUM})",re.I)
    tmp={}
    for line in path.read_text().splitlines():
        m=pat.match(line)
        if m:
            key=(int(m.group(2)),int(m.group(3)))
            tmp.setdefault(key,{})[m.group(1).lower()]=norm(m.group(4))
    return {k:(v.get('alpb','0'),v.get('xfac','0')) for k,v in tmp.items() if 'alpb' in v and 'xfac' in v}

def write_csv(path:Path,sources:list[str],header,rows):
    with path.open('w',newline='',encoding='utf-8') as f:
        f.write('# PROVENANCE: extracted from OpenMOPAC (Molecular Orbital PACkage)\n')
        f.write(f'# source repo: {REPO}  tag: {TAG}\n')
        for s in sources: f.write(f'# source file: {s}\n')
        f.write('# source license: Apache-2.0\n')
        f.write('# extraction script: tools/extract_xndo_params.py\n')
        w=csv.writer(f,lineterminator='\n'); w.writerow(header); w.writerows(rows)
    print(f'wrote {path} ({len(rows)} rows)')

def main():
    ap=argparse.ArgumentParser();
    ap.add_argument('--mndo',type=Path,required=True); ap.add_argument('--mndod',type=Path,required=True)
    ap.add_argument('--mndo-pairs',type=Path,required=True); ap.add_argument('--mndod-pairs',type=Path,required=True)
    ap.add_argument('--indo',type=Path,required=True); args=ap.parse_args()
    out=Path(__file__).resolve().parents[1]/'src'/'data'; out.mkdir(exist_ok=True)
    cols=['z','sym','uss','upp','udd','zs','zp','zd','betas','betap','betad','gss','gsp','gpp','gp2','hsp','zsn','zpn','zdn','f0sd','g2sd','alp','poc','g1_k','g1_l','g1_m','g2_k','g2_l','g2_m','g3_k','g3_l','g3_m','g4_k','g4_l','g4_m']
    def nddo(path:Path,kind:str):
        d=parse_data(path)
        if kind=='mndo':
            mp={'uss':'ussm','upp':'uppm','udd':'uddm','zs':'zsm','zp':'zpm','zd':'zdm','betas':'betasm','betap':'betapm','betad':'betadm','gss':'gssm','gsp':'gspm','gpp':'gppm','gp2':'gp2m','hsp':'hspm','zsn':'zsnm','zpn':'zpnm','zdn':'zdnm','f0sd':'f0sdm','g2sd':'g2sdm','alp':'alpm','poc':'pocm'}
            gs=('guesm1','guesm2','guesm3')
        else:
            mp={'uss':'ussd','upp':'uppd','udd':'uddd','zs':'zsd','zp':'zpd','zd':'zdd','betas':'betasd','betap':'betapd','betad':'betadd','gss':'gssd','gsp':'gspd','gpp':'gppd','gp2':'gp2d','hsp':'hspd','zsn':'zsnd','zpn':'zpnd','zdn':'zdnd','f0sd':None,'g2sd':None,'alp':'alpd','poc':'poc_d'}
            gs=(None,None,None)
        def get(name,z,k=None):
            if name is None:return '0'
            return d.get((name,z,k),'0')
        zs=sorted({z for (name,z,k) in d if name in {v for v in mp.values() if v}})
        rows=[]
        for z in zs:
            vals=[z,SYMBOLS[z-1] if 1<=z<=len(SYMBOLS) else f'Z{z}']
            for c in cols[2:23]: vals.append(get(mp[c],z))
            for g in range(1,5):
                for gn in gs: vals.append(get(gn,z,g))
            rows.append(vals)
        return rows
    mndo_rows=nddo(args.mndo,'mndo'); mndod_rows=nddo(args.mndod,'mndod')
    assert abs(float(mndo_rows[0][2])-(-11.906276))<1e-10
    write_csv(out/'mndo_parameters.csv',['src/models/parameters_for_mndo_C.F90'],cols,mndo_rows)
    write_csv(out/'mndod_parameters.csv',['src/models/parameters_for_mndod_C.F90'],cols,mndod_rows)
    for name,path in [('mndo',args.mndo_pairs),('mndod',args.mndod_pairs)]:
        pairs=parse_pairs(path); rows=[[a,b,x,y] for (a,b),(x,y) in sorted(pairs.items())]
        write_csv(out/f'{name}_pair_parameters.csv',[f'src/models/alpb_and_xfac_{name}.F90'],['zi','zj','alpb','xfac'],rows)
    # ZINDO/S table: preserve all 24 Zerner FG parameters plus double-zeta d data.
    z=parse_data(args.indo)
    zh=['z','sym','isok','norb','zcore','zeta_sp','zeta_d1','zeta_d2','zeta_w1','zeta_w2','beta_s','beta_p','beta_d']+[f'fg{i}' for i in range(1,25)]
    rows=[]
    for zz in range(1,81):
        isok=int(float(z.get(('isoki',zz,None),'0')))
        norb=int(float(z.get(('nbfai',zz,None),'0')))
        if not isok: continue
        row=[zz,SYMBOLS[zz-1],isok,norb,z.get(('zcoreai',zz,None),'0'),z.get(('zetai',zz,None),'0'),z.get(('zetadi',1,zz),'0'),z.get(('zetadi',2,zz),'0'),z.get(('zetawti',1,zz),'0'),z.get(('zetawti',2,zz),'0'),z.get(('betaai',1,zz),'0'),z.get(('betaai',2,zz),'0'),z.get(('betaai',3,zz),'0')]
        row += [z.get(('fgi',i,zz),'0') for i in range(1,25)]
        rows.append(row)
    write_csv(out/'zindo_s_parameters.csv',['src/models/parameters_for_INDO_C.F90'],zh,rows)
    print('ZINDO/S supported elements:',len(rows),[r[1] for r in rows])

if __name__=='__main__': main()
