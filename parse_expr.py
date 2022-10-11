#!/usr/bin/env python3
import os
import csv

def Average(lst):
    return sum(lst) / len(lst)

diseases=['ibd','pso','ra','sle']
dir=os.environ['TABT_DIR']
os.chdir(dir+'/data/expression')

for dis in diseases:
    with open(dis+'_expr_files.tsv') as f:
        efns=f.readlines()
        expr_files=[efn.rstrip() for efn in efns]

    
    tes=dict()
    gene_count=dict()
    for efn in expr_files:
        tes[efn]=dict()
        with open(dis+'/'+efn+'.top.table.tsv') as f:
            reader = csv.reader(f, delimiter='\t')
            next(reader)
            limma=list(reader)
        for li in limma:
            t=li[3]
            g=li[6]
            if (g!='') and (g not in tes[efn]): 
                tes[efn][g]=t
        for g in tes[efn].keys():
            if g not in gene_count:
                gene_count[g]=1
            else:
                gene_count[g]+=1
    with open(dis+'_consensus.tsv','w') as f:    
        print('Gene',end='\t',file=f)
        for efn in expr_files: print(efn,end='\t',file=f)
        print('Average',file=f)
        for g in gene_count.keys():
            if gene_count[g] > len(efns)/2:
                print(g,end='\t',file=f)
                avel=list()
                for efn in expr_files:
                    try:
                        print(tes[efn][g],end='\t',file=f)
                        avel.append(float(tes[efn][g]))
                    except:
                        print('NA',end='\t',file=f)
                print(Average(avel),file=f)
