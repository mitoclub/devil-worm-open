## by hand: 
## rename unitig_62_quiver_annotated_Trimmed_CURATED_ANNOTATIONS_..tsv into unitig_62_quiver_annotated_Trimmed_CURATED_ANNOTATIONS_Halicephalobus mephisto.tsv
## pachys_mito_Bracht_Curated_synthetic DNA construct_..tsv into pachys_mito_Bracht_Curated_synthetic DNA construct_Diploscapter pachys.tsv

rm(list=ls(all=TRUE))

CU <- read.table("../../Body/2Derived/CodonUsageDerived.csv", sep = ';', head = TRUE)
table(CU$OrganismName)
length(unique(CU$OrganismName)) # 37
table(CU$Strand)  # 
table(CU$GenName) # should be 13
#atp6       ATP6       ATP8        cob       cox1       Cox1       COX1       cox2       COX2       cox3       COX3       CYTB  Dpa-atp-6 
#1         36          7          1          1          1         33          1         34          1         34         34          1 
#Dpa-ctb-1  Dpa-ctc-3 Dpa-ndfl-4 Dpa-nduo-1 Dpa-nduo-2 Dpa-nduo-5 Dpa-nduo-6       nad1       NAD1       nad2       NAD2       nad3       NAD3 
#1          1          1          1          1          1          1          1          2          1          2          1          2 
#nad4       NAD4      nad4L      NAD4L       nad5       NAD5       nad6       NAD6        ND1        ND2        ND3        ND4       ND4L 
#1          2          1          2          1          2          1          2         32         32         36         33         32 
#ND5        ND6       None 
#32         32         12

#### -1 strand (who they are?) if we are sure => we just reverse compliment in the next step? For now we delete them
CU = CU[CU$Strand == 1,]
table(CU$OrganismName)

#### aggregate all data by species
names(CU)
CUAggByGenes = aggregate(CU[8:71], by = list(CU$GenName), FUN = sum) # as soon as gene names are good => check it better
CUAggBySpecies = aggregate(CU[8:71], by = list(CU$OrganismName), FUN = sum) # as soon as gene names are good => check it better
names(CUAggBySpecies)[1]='OrganismName'

### derive skew towards T for fourfold degenerate sites (expectation is 25%) (assume invertebrate genetic code, will remove homo sapiens later):

CUAggBySpecies$SerAG.T = CUAggBySpecies$AGT/(CUAggBySpecies$AGA+CUAggBySpecies$AGT+CUAggBySpecies$AGG+CUAggBySpecies$AGC) # AGX = Ser.AG
CUAggBySpecies$SerTC.T = CUAggBySpecies$TCT/(CUAggBySpecies$TCA+CUAggBySpecies$TCT+CUAggBySpecies$TCG+CUAggBySpecies$TCC) # AGX = Ser.TC
CUAggBySpecies$Leu.T = CUAggBySpecies$CTT/(CUAggBySpecies$CTA+CUAggBySpecies$CTT+CUAggBySpecies$CTG+CUAggBySpecies$CTC) # CTX = Leu
CUAggBySpecies$Pro.T = CUAggBySpecies$CCT/(CUAggBySpecies$CCA+CUAggBySpecies$CCT+CUAggBySpecies$CCG+CUAggBySpecies$CCC) # CCX = Pro
CUAggBySpecies$Arg.T = CUAggBySpecies$CGT/(CUAggBySpecies$CGA+CUAggBySpecies$CGT+CUAggBySpecies$CGG+CUAggBySpecies$CGC) # CGX == Arg (R)
CUAggBySpecies$Thr.T = CUAggBySpecies$ACT/(CUAggBySpecies$ACA+CUAggBySpecies$ACT+CUAggBySpecies$ACG+CUAggBySpecies$ACC) # ACX == Thr
CUAggBySpecies$Val.T = CUAggBySpecies$GTT/(CUAggBySpecies$GTA+CUAggBySpecies$GTT+CUAggBySpecies$GTG+CUAggBySpecies$GTC) # GTX == Val
CUAggBySpecies$Ala.T = CUAggBySpecies$GCT/(CUAggBySpecies$GCA+CUAggBySpecies$GCT+CUAggBySpecies$GCG+CUAggBySpecies$GCC) # GCX == Ala
CUAggBySpecies$Gly.T = CUAggBySpecies$GGT/(CUAggBySpecies$GGA+CUAggBySpecies$GGT+CUAggBySpecies$GGG+CUAggBySpecies$GGC) # GGX == Gly (G)

### derive mean skew and sort
names(CUAggBySpecies)
T = as.data.frame(CUAggBySpecies[,66:74])
CUAggBySpecies$SkewTowardsT4Dsites = apply(T, 1, mean)
CUAggBySpecies = CUAggBySpecies[order(CUAggBySpecies$SkewTowardsT4Dsites),]
summary(CUAggBySpecies$SkewTowardsT4Dsites)
Champions = CUAggBySpecies[CUAggBySpecies$SkewTowardsT4Dsites > 0.75,]$OrganismName
Champions

### overlap these results with the tree!!! 

### what about Amino Acids? 


