import sys,pandas

INFO_FILE=sys.argv[1]

tf2fam,mot2fam,mot2tf={},{},{}
with open(INFO_FILE) as f:
    for line in f:
        spl=line.strip().split('\t')
        tf2fam[spl[0]]=spl[1]
        if spl[2]=='Y':
            for m in spl[3].rsplit(','):
                if m in mot2fam:
                    mot2fam[m]+=[spl[1]]
                    mot2tf[m]+=[spl[0]]
                else:
                    mot2fam[m]=[spl[1]]
                    mot2tf[m]=[spl[0]]
for m in mot2fam:
    mot2fam[m]=list(set(mot2fam[m]))
    mot2tf[m]=list(set(mot2tf[m]))
    


def enricher_parser(infile,outfile,extended):
    fin=open(infile,'r')
    df=[]
    for line in fin:
        if line.startswith('#') or line.startswith('set_id'): # skip the comments and the header
            pass
        else:
            spl=line.strip().split("\t")
            target_genes = spl[9].rsplit(",")
            if extended == 'TF_motifs':
                if spl[0] in mot2tf[spl[1]]:
                    for g in target_genes:
                        df.append([spl[0],g])
            elif extended == 'TF-F_motifs':
                if spl[0] in tf2fam and tf2fam[spl[0]] in mot2fam[spl[1]]:
                    if tf2fam[spl[0]] != 'Unknown':
                        for g in target_genes:
                            df.append([spl[0],g])  
                    else:
                        if spl[0] in mot2tf[spl[1]]: #do not extend to fams when the fam is "Unknown"
                            for g in target_genes:
                                df.append([spl[0],g]) 


    df=pandas.DataFrame(df)
    df=df.drop_duplicates().reset_index(drop=True)
    df.to_csv(outfile,sep='\t',index=None,header=None)
       
inf =sys.argv[2]
outf = sys.argv[3]
ext = sys.argv[4]
enricher_parser(inf,outf,ext)