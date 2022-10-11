#!/usr/bin/env python3
import os
import csv
diseases=['ibd','pso','ra','sle']

dir=os.environ['TABT_DIR']
os.chdir(dir+'/data/gwas')

for dis in diseases:
    with open(dis+'.tsv') as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader)
        gwas = list(reader)
        f.close() 
    
    cgwas=dict()
    for li in gwas:
        pval=float(li[27])
        if (pval <= 1e-20):
            mgenes=li[14].replace(' - ',', ').replace(' x ',', ').split(',')
            for ge in mgenes:
                ge=ge.strip()
                if (ge=='NR'): continue
                if (ge==''): continue
                if (ge.lower()=='intergenic'): continue
                if (ge.lower()=='no mapped genes'): continue
                if (ge in cgwas):
                    cgwas[ge]+=1  
                else: 
                    cgwas[ge]=1 
                    
    with open(dis+'_genes.tsv','w') as f:
        for ge in cgwas.keys():
            print(ge,cgwas[ge],sep='\t',file=f)
        f.close() 
