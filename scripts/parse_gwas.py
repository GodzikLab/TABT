#!/usr/bin/env python3
import os
import csv
diseases=['ibd','pso','ra','sle']

dir=os.environ['TABT_DIR']
os.chdir(dir+'/data/gwas')

for dis in diseases:
    print('Reading GWAS data for '+dis+'...')
    with open(dis+'.tsv') as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader)
        gwas = list(reader)
        f.close() 
    print('Parsing GWAS data for '+dis+'...')
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
    print('Saving list of genes with mutations associated with '+dis)              
    with open(dis+'_genes.tsv','w') as f:
        for ge in cgwas.keys():
            print(ge,cgwas[ge],sep='\t',file=f)
        f.close() 
print('Done')
