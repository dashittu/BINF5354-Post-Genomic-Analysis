# parsevcf.py
# Annotates variants from CSV files with gene names, region types, and mutation types
# Original code by Dayo Shittu

from math import fabs
import os
import pandas as pd
import gc

AADict = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"*", "TAG":"*",
    "TGT":"C", "TGC":"C", "TGA":"*", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",}
def annotation(unify,MAX_DIST = 5000,n=3):# unify is the csv file you generated in homework 2
    dictionary= {} 
    position = 0
    chrom=''
    prechrom=0
    
    for l in range(unify.shape[0]): #line in vcf:
        chrom=unify["chrom"][l]
        position = unify["left"][l]
        ref=unify["ref_seq"][l]
        if unify["ref_seq"][l]==unify["var_seq1"][l]:
            alt=unify["var_seq2"][l] 
        else:
            alt=unify["var_seq1"][l]
        gtfCounter = -1 #point to each line in gtf file
        genestart = 0
        genestop = 0
        boo = 0
        where = []
        genename = []
        change_type=[]
        counter=[]
        strand_gene=[]
        if chrom != prechrom:
            # You need to change this and the following directory
            f = open("/Users/dayoshittu/Documents/Fall2022/BINF5354/ref_genome//"+chrom+'_hg38refFlat_2020.txt','r')
            flat = f.readlines()
            f.close()
            seqfile=open("/Users/dayoshittu/Documents/Fall2022/BINF5354/ref_genome/genome_"+chrom+".fa",'r')
            next(seqfile)
            sequence=seqfile.read().replace('\n','')
            seqfile.close()
        prechrom=chrom
        
        while boo == 0 and gtfCounter < len(flat):
            if position > genestart and position <= genestop and gtfCounter < (len(flat)-1): #GSV within gene
                counter.append(gtfCounter)
                strand_gene.append(flat[gtfCounter].split()[3])
                gtfCounter += 1
                pgstop = genestop
                genestart = int(flat[gtfCounter].split()[4]) #transcript start
                genestop = int(flat[gtfCounter].split()[5]) #transcript end
            elif position > genestop and gtfCounter < (len(flat)-1):  #Moves to the next gene
                gtfCounter += 1
                pgstop = genestop #prev gene stop
                genestart = int(flat[gtfCounter].split()[4]) #transcript start
                genestop = int(flat[gtfCounter].split()[5]) #transcript end
            elif position < genestart and position > pgstop:
                boo = 1                       
            else:
                boo = 1
        if len(strand_gene)==0:
            strand_gene.append('')
        strand = flat[gtfCounter].split()[3] #ref strand direction
        if len(counter)==0:
            if gtfCounter == len(flat)-1 and position > genestop:
                if (position - genestop) > MAX_DIST: #if variant position is far from the end of transcript
                    genename.append('NoName')
                    where.append('Not close to gene')
                    change_type.append('')
                else:
                    genename.append(flat[gtfCounter].split()[0]) #else whether close to 3' or 5'
                    if strand == '-':
                        where.append('5\' untranscribed')
                        change_type.append('')
                    else:
                        where.append('3\' untranscribed')
                        change_type.append('')                        
            elif position <= genestart and position > pgstop:
                # check if closer to start or end of the one before or current
                diff1 = fabs(position-genestart) #current
                diff2 = fabs(position-pgstop) #previous
                if diff1 <= diff2 and diff1 <= MAX_DIST: #if closer to current "new " gene and within cutoff distance
                    genename.append(flat[gtfCounter].split()[0])
                    if strand == '-':
                        where.append('3\' untranscribed')
                        change_type.append('')
                    else:
                        where.append('5\' untranscribed')
                        change_type.append('')                        
                elif diff2 <= MAX_DIST: #if closer to prev gene and within cutoff distance
                    genename.append(flat[gtfCounter-1].split()[0])
                    strand = flat[gtfCounter-1].split()[3]
                    if strand == '-':
                        where.append('5\' untranscribed')
                        change_type.append('')
                    else:
                        where.append('3\' untranscribed')
                        change_type.append('')                       
                else:
                    genename.append('NoName')
                    where.append('Not close to gene')
                    change_type.append('')
            else:
                #continue
                print("next variant")
        else:

            for i in range(len(counter)):
                genename.append(flat[counter[i]].split()[0])
                cd_start=int(flat[counter[i]].split()[6]) #CDS start
                cd_end=int(flat[counter[i]].split()[7]) #CDS end
                ex_starts = []
                ex_ends=[]
                x=0
                tmp = 'Unknown'
                ex_starts = [int(i) for i in flat[counter[i]].split()[-2].split(',')[0:-1]]
                ex_ends = [int(i) for i in flat[counter[i]].split()[-1].split(',')[0:-1]]
    
                while tmp == 'Unknown' and x < len(ex_starts):
                    ######
                    if position > ex_starts[x] and position <= ex_ends[x]:  # if position in exon
                        tmp = 'Exon'
                    elif position <= ex_starts[x]:
                        tmp = 'Intron'
                    else:
                        x += 1
                if tmp == 'Intron' or tmp=='Unknown':
                    where.append('Intron')
                    change_type.append('')
                #see if position is in translated region:
                elif tmp == 'Exon':
                    if position >  cd_start and position <= cd_end:
                        tmp = 'Translated'
                    else: #not translated
                        tmp = 'Not translated'
                    if tmp == 'Not translated':
                        if cd_start == cd_end: #ncRNAs
                            where.append('Not translated, ncRNA')
                            change_type.append('')
                        else: #utrs
                            if position  <= cd_start:
                                if strand_gene[i] == '-':
                                    where.append('Not translated,3\' utr')
                                    change_type.append('')
                                else:
                                    where.append('Not translated,5\' utr')
                                    change_type.append('')
                            else:
                                if strand_gene[i] == '-':
                                    where.append('Not translated,5\' utr')
                                    change_type.append('')
                                else:
                                    where.append('Not translated,3\' utr')
                                    change_type.append('')
                    elif tmp == 'Translated':
                        where.append('CDS')
                        ref_codon=''
                        var_codon=''
                        m=0
                        toRev=0
                        if strand_gene[i]=="+":
                            ex=0
                            while ex != (len(ex_ends)-1) and position > ex_ends[ex]:
                                if cd_start > ex_ends[ex]:
                                    ex+=1
                                elif cd_start > ex_starts[ex] and cd_start < ex_ends[ex] and position > ex_ends[ex]:
                                    m = ((ex_ends[ex] - cd_start )+m)%3
                                    ex+=1
                                else:
                                    m = ((ex_ends[ex]-ex_starts[ex])+m)%3
                                    ex+=1
                            if cd_start > ex_starts[ex]:
                                m=(position-cd_start-1+m)%3
                            else:
                                m=(position-ex_starts[ex]-1+m)%3
            
                            if m==0 :
                                ref_codon=ref+sequence[position]+sequence[position+1]
                                var_codon=alt+sequence[position]+sequence[position+1]
                            elif m==1:
                                ref_codon = sequence[position-2]+ref+ sequence[position]
                                var_codon = sequence[position - 2]+alt+sequence[position]
                            else:
                                ref_codon = sequence[position-3]+sequence[position-2]+ref
                                var_codon = sequence[position-3]+sequence[position-2]+alt                                                             
                        else:
                            ex = len(ex_starts)-1
                            while ex != 0 and position < int(ex_starts[ex]):
                                #######
                                if cd_end < ex_starts[ex]:
                                    ex -= 1
                                elif cd_end < ex_ends[ex] and cd_end > ex_starts[ex] and position < ex_starts[ex]:
                                    m = ((cd_end - ex_starts[ex]) + m) % 3
                                    ex -= 1
                                else:
                                    m = ((ex_ends[ex] - ex_starts[ex]) + m) % 3
                                    ex -= 1
                            
                            if cd_end > ex_ends[ex]:
                                m=(ex_ends[ex] - position + m)%3
                            else:
                                m=(cd_end - position + m)%3

                            toRev=1
                            if m==0 :
                                ref_codon=ref+sequence[position-2]+sequence[position-3]
                                var_codon=alt+sequence[position-2]+sequence[position-3]
                            elif m==1:
                                ref_codon = sequence[position]+ref+ sequence[position-2]
                                var_codon = sequence[position]+alt+sequence[position-2]
                            else:
                                ref_codon = sequence[position+1]+sequence[position]+ref
                                var_codon = sequence[position+1]+sequence[position]+alt
                            
                        if toRev==1:
                            intab = "ATCG"
                            outtab = "TAGC"
                            trantab = str.maketrans(intab, outtab)
                            var_codon = var_codon.translate(trantab)
                            ref_codon = ref_codon.translate(trantab)
                            
                        var_codon = var_codon[0:3]
                        ref_codon = ref_codon[0:3]
                        ref_aa=''
                        var_aa=''
                        if var_codon in AADict:
                            var_aa = AADict[var_codon] #translate
                        if ref_codon in AADict:
                            ref_aa = AADict[ref_codon]
                        if var_aa == ref_aa:
                            change_type.append('Synonymous')
                        else:
                            change_type.append('Non-synonymous')
                    else:  #To make sure to fill change type in on the ignored ones.
                        where.append('Ignore')
                        change_type.append('')
                else:
                    where.append('Ignore')
                    change_type.append('')
        dictionary[l]=[chrom,position,strand,genename,where,change_type] 
    del flat,seqfile
    gc.collect()
    return dictionary 

def combine(geneinfo,csvfile):
    genename=[]
    where=[]
    change_type=[]
    for value in geneinfo.values():
        '''
        Delete this comment and fill in the block. Append values for each variant to corresponding list.        
        '''
        ####### MY CODE GOES HERE ###########
        genename.append(value[3][0])
        where.append(value[4][0])
        change_type.append(value[5][0])
        
    
# The following codes insert three columns to csv file and save the dataframe to the same file(overwrite)
# If you want to save to different directory, you may add one more parameter to the function to specify output path         
    csvf=pd.read_csv(csvfile)    
    csvf.insert(csvf.shape[1],"genename",genename)
    csvf.insert(csvf.shape[1],"where",where)
    csvf.insert(csvf.shape[1],"change_type1",change_type)  
    csvf.to_csv(csvfile,index=False)
    return


if __name__=='__main__':
    os.chdir("/Users/dayoshittu/Documents/Fall2022/BINF5354/Project/merged")
    csv_files = os.listdir()
    
    for file in csv_files:
        unify = annotation(pd.read_csv(file), 5000, 3)
        combine(unify, file)
    
        
        
        
        
        
        
        
        
