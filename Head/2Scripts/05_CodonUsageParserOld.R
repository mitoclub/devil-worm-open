## by hand: 
## rename unitig_62_quiver_annotated_Trimmed_CURATED_ANNOTATIONS_..tsv into unitig_62_quiver_annotated_Trimmed_CURATED_ANNOTATIONS_Halicephalobus mephisto.tsv
## pachys_mito_Bracht_Curated_synthetic DNA construct_..tsv into pachys_mito_Bracht_Curated_synthetic DNA construct_Diploscapter pachys.tsv

rm(list=ls(all=TRUE))

Final=data.frame()
Dirs <- dir("../../Body/2Derived/CodonUsageOld/") # /home/popadin/devil-worm/Body/2Derived
for (i in 1 : length(Dirs)) 
{ # i = 1
  Path=paste('../../Body/2Derived/CodonUsageOld/',Dirs[i],sep='')
  if (file.exists(Path))
  {
    Species = gsub('(.*)_','',Path); Species = gsub('.tsv','',Species); Species = gsub(' ','_',Species);    
    if (Species != '.')
     {
     sp = read.table(Path, head = TRUE)
     
     #### RENAME
     VecOfOldNames = colnames(sp)
     VecOfNewNames=c()
     for (j in 1:length(VecOfOldNames))
      { # j = 1
        temp = gsub('(.*)\\.','',VecOfOldNames[j])
        temp = toupper(temp)
        VecOfNewNames = c(VecOfNewNames,temp)
      }
     names(sp)=VecOfNewNames
     
     #### COUNT TOTAL CODON USAGE AND SAVE TWO-COLUMN FILE:
    
      sp$Total = 0
      for (k in 1:nrow(sp)) { sp$Total[k]=  sum(sp[k,]) }
      sp$Codons = row.names(sp)
      sp=sp[colnames(sp) %in% c('Total','Codons')]
      names(sp)[1]=Species;
    
     #### MERGE WITH FINAL
     if (i == 1) {Final = sp}
     if (i > 1) {Final = merge(Final,sp, by = 'Codons')} 
    }
  }
}  

write.table(Final,"../../Body/3Results/01.CodonUsageParserOld.TotalCodonTable.txt")

########## ANALYSIS OF THE TOTAL CODON USAGE
###### ONLY FOUR FOLD DEGENER SITES FROM GENETIC CODE 5

VecOfFourFoldDegSites = c(
'TCA','TCT','TCG','TCC',
'CTA','CTT','CTG','CTC',
'CCA','CCT','CCG','CCC',
'CGA','CGT','CGG','CGC',
'ACA','ACT','ACG','ACC',
'AGA','AGT','AGG','AGC',
'GTA','GTT','GTG','GTC',
'GCA','GCT','GCG','GCC',
'GGA','GGT','GGG','GGC')

FourFold = Final[Final$Codons %in% VecOfFourFoldDegSites,]
FourFold = FourFold[order(FourFold$Codons),]
# even by eye we can see an excess of T and A

LastA = c('TCA','CTA','CCA','CGA','ACA','AGA','GTA','GCA','GGA')
LastT = c('TCT','CTT','CCT','CGT','ACT','AGT','GTT','GCT','GGT')
LastG = c('TCG','CTG','CCG','CGG','ACG','AGG','GTG','GCG','GGG')
LastC = c('TCC','CTC','CCC','CGC','ACC','AGC','GTC','GCC','GGC')

FourFoldLastA = Final[Final$Codons %in% LastA,]; 
FourFoldLastA = FourFoldLastA[!colnames(FourFoldLastA) %in% 'Codons']; 
A=data.frame(apply(FourFoldLastA,2,sum)); names(A)=c('A')

FourFoldLastT = Final[Final$Codons %in% LastT,]; 
FourFoldLastT = FourFoldLastT[!colnames(FourFoldLastT) %in% 'Codons']; 
T=data.frame(apply(FourFoldLastT,2,sum)); names(T)=c('T')

FourFoldLastG = Final[Final$Codons %in% LastG,]; 
FourFoldLastG = FourFoldLastG[!colnames(FourFoldLastG) %in% 'Codons']; 
G=data.frame(apply(FourFoldLastG,2,sum)); names(G)=c('G')

FourFoldLastC = Final[Final$Codons %in% LastC,]; 
FourFoldLastC = FourFoldLastC[!colnames(FourFoldLastC) %in% 'Codons']; 
C=data.frame(apply(FourFoldLastC,2,sum)); names(C)=c('C')

ATGC = cbind(A,T,G,C)
ATGC$FrA = ATGC$A/(ATGC$A+ATGC$T+ATGC$G+ATGC$C)
ATGC$FrT = ATGC$T/(ATGC$A+ATGC$T+ATGC$G+ATGC$C)
ATGC$FrG = ATGC$G/(ATGC$A+ATGC$T+ATGC$G+ATGC$C)
ATGC$FrC = ATGC$C/(ATGC$A+ATGC$T+ATGC$G+ATGC$C)
ATGC$FrAT = (ATGC$A+ATGC$T)/(ATGC$A+ATGC$T+ATGC$G+ATGC$C)
ATGC$FrGC = (ATGC$G+ATGC$C)/(ATGC$A+ATGC$T+ATGC$G+ATGC$C)

ATGC = ATGC[order(ATGC$FrT),]

write.table(ATGC,"../../Body/3Results/01.CodonUsageParserOld.ATGCFreqAtFourFoldDegSites.txt")

### an excess of T in mephisto can be very low A>G on opposite chain, or high C>T on light strand
### all low-Ne species becomes AT rich (endosymbisys => no reparatio)


### AMINOACIDS:

#TCX Ser1
#CTX Leu
#CCX Pro
#CGX Arg
#ACX Thr
#AGX Ser2
#GTX Val
#GCX Ala
#GGX Gly

###### GA SKEW

