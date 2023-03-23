#!/usr/bin/env python3
import os
import csv
import urllib.request
import gzip

def Average(lst):
    return sum(lst) / len(lst)

diseases=['ibd','pso','ra','sle']

dir=os.environ['TABT_DIR']

print('Preparing protein-protein interactions data...')
os.chdir(dir+'/data/interactions')
print("Downloading STRING datasets...")

urllib.request.urlretrieve("https://stringdb-static.org/download/protein.aliases.v11.5/9606.protein.aliases.v11.5.txt.gz","aliases.gz")
urllib.request.urlretrieve("https://stringdb-static.org/download/protein.links.v11.5/9606.protein.links.v11.5.txt.gz","links.gz")


print("Reading protein-protein interactions data...")
with gzip.open('links.gz','rt') as f:
    reader = csv.reader(f, delimiter=" ")
    sinter = list(reader)

with gzip.open('aliases.gz','rt') as f:
    reader = csv.reader(f, delimiter="\t")
    salias = list(reader)

print("Translating gene names...")
galias=dict()
cint=dict()

for li in salias:
    if (li[2]=='BioMart_HUGO'): galias[li[0]]=li[1]

for li in salias:
    if ((li[2]=='Ensembl_UniProt_GN') and (li[0] not in galias)): galias[li[0]]=li[1]
    


ints=dict()
for li in sinter:
    if ((li[0] in galias) and (li[1] in galias)):
        ge0=galias[li[0]]
        ge1=galias[li[1]]
        if (ge1 < ge0):  gt=ge1; ge1=ge0; ge0=gt
        score=li[2]
        if (ge0 not in ints): ints[ge0]=dict()
        if (ge1 not in ints[ge0]):  ints[ge0][ge1]=score

f = open("interactions.tsv", "w")
for ge0 in ints:
    print('>',ge0,sep='',file=f)
    for ge1 in ints[ge0]:
        print(ge1,file=f)
f.close()
print("Done")

print('Preparing differential expression data...')
os.chdir(dir+'/data/expression')

for dis in diseases:
    print('Reading differential expression data for '+dis+'...')
    with open(dis+'_expr_files.tsv') as f:
        efns=f.readlines()
        expr_files=[efn.rstrip() for efn in efns]

    print('Preparing differential expression consensus list for '+dis+'...')
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
    print('Saving differential expression consensus list for '+dis+'...')
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
print('Done')

print('Preparing GWAS data...')
os.chdir(dir+'/data/gwas')

for dis in diseases:
    print('Reading GWAS data for '+dis+'...')
    with open(dis+'.tsv') as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader)
        gwas = list(reader)
        f.close() 
    print('Filtering GWAS data for '+dis+'...')
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