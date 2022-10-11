#!/usr/bin/env python3
import os
import csv

dir=os.environ['TABT_DIR']
os.chdir(dir+'/data/interactions')

with open('9606.protein.links.v11.5.txt') as f:
    reader = csv.reader(f, delimiter=" ")
    sinter = list(reader)

with open('9606.protein.aliases.v11.5.txt') as f:
    reader = csv.reader(f, delimiter="\t")
    salias = list(reader)

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

        