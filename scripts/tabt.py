#!/usr/bin/env python3

import os
import csv

diseases=['ibd','pso','ra','sle']

dir=os.environ['TABT_DIR']
os.chdir(dir)

def Tabt(gs,inv,interd,interc,invin,prop,tabt):
    for g in gs:
        invin[g]=0
        for ig in interd[g]:
            if ig in inv: 
                invin[g]+=1
        prop[g]=invin[g]/interc[g]
        
    propave=sum(prop.values())/len(prop)
        
    for g in gs:
        if (interc[g]>10):
            tabt[g]=(prop[g]-propave)/propave
        else:tabt[g]=0
    return 1


with open('data/location/uniprot.tsv') as f:
    reader = csv.reader(f, delimiter='\t')
    next(reader)
    unip = list(reader)

local=dict()
for uli in unip:
    g=uli[6]
    l=uli[5].lower()
    if (('membrane' in l) or ('secreted' in l)): local[g]=1
    else: local[g]=0
    
with open('data/interactions/interactions.tsv') as f:
    reader = csv.reader(f, delimiter='\t')
    interl = list(reader)

interd=dict()
interc=dict()

for ili in interl:
    g=ili[0]
    if ('>' in g):
        ge0=g.replace('>','')
        if (ge0 not in interd): interd[ge0]=list()
    else: 
        ge1=g
        if (ge1 not in interd): interd[ge1]=list()
        interd[ge0].append(ge1)
        interd[ge1].append(ge0)

for g in interd: interc[g]=len(interd[g])

gset=(local.keys() & interc.keys())

propll=dict()

with open('data/tissues/E-MTAB-513-query-results.tpms.tsv') as f:
    reader = csv.reader(filter(lambda row: row[0]!='#',f), delimiter='\t')
    next(reader)
    tiss=list(reader)
        
for tli in tiss:
    if (tli[1] in gset):
        tlif=list(map(lambda x: float(x or 0), tli[2:]))
        propll[tli[1]]=(tlif[7]+tlif[10])/sum(tlif)

propllave=sum(propll.values())/len(gset)

tabtll=dict()

for g in gset:
    if (g not in propll): propll[g]=0
    tabtll[g]=(propll[g]-propllave)/propllave

tabt=dict()
for dis in diseases:
    tabt[dis]=dict()
    with open('data/expression/'+dis+'_consensus.tsv') as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader)
        consl=list(reader)
    invex=dict()
    for cli in consl:
        g=cli[0]
        ave=cli[-1]
        if (abs(float(ave)) > 5): invex[g]=1
        
    with open('data/gwas/'+dis+'_genes.tsv') as f:
        reader = csv.reader(f, delimiter='\t')
        gwas=list(reader)
    invgs=dict()
    for gli in gwas:
        invgs[gli[0]]=1
    

    invinex=dict(); propex=dict(); tabtex=dict();    
    invings=dict(); propgs=dict(); tabtgs=dict();
    
    Tabt(gset, invex, interd, interc, invinex, propex, tabtex)
    Tabt(gset, invgs, interd, interc, invings, propgs, tabtgs)
    

    for g in gset:
        tabt[dis][g]=max(local[g]*(tabtll.get(g)+tabtgs.get(g)*0+tabtex.get(g)),0)
    
    with open('results/'+dis+'_tabt.tsv','w') as f: 
        for g in gset:
            print(g,interc[g],
                  invinex[g],propex[g],tabtex[g],
                  invings[g],propgs[g],tabtgs[g],
                             propll[g],tabtll[g],
                              local[g],tabt[dis][g],sep='\t',file=f)


tpercselfave=0
with open('results/evaluation.tsv','w') as f:
    for dis in diseases:
        tlen=len(tabt[dis])
        tabtperc={k:i/tlen for i,k in 
                  enumerate(sorted(tabt[dis],key=tabt[dis].get,reverse=True))}
        for distst in diseases:
            with open('data/abs/'+distst+'_abs.tsv') as g:
                reader = csv.reader(g, delimiter='\t')
                tabs=list(reader)
                tpercave=0
            for tab in tabs: 
                print(dis,distst,tab[0],tabt[dis][tab[0]],
                      tabtperc[tab[0]],sep='\t',file=f)
                tpercave+=tabtperc[tab[0]]
            tpercave/=len(tabs)
            if (dis==distst): tpercselfave+=tpercave
            print('average:',dis,distst,tpercave,sep='\t',file=f)
    tpercselfave/=len(diseases)
    
    print('Average:',tpercselfave,sep='\t',file=f)
    
        
        