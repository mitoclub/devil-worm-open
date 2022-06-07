WholeGenomeContent <- read.table('/home/emulciber/MitoClub/devil-worm/Body/3Results/10_GenomeStatistics.csv', head = TRUE, sep = ';')
FourFoldContent <- read.table('/home/emulciber/MitoClub/devil-worm/Body/3Results/09_FourFoldContent.csv', head = TRUE, sep = ';')

WholeGenomeContent[c(1,3:7)] <- list(NULL)
FourFoldContent[c(2:65)] <- list(NULL)

rownames(WholeGenomeContent) <- WholeGenomeContent$Species
rownames(FourFoldContent) <- FourFoldContent$OrganismName

Species <- which(colnames(WholeGenomeContent) == c("Species", "Group"))
SpeciesF <- which(colnames(FourFoldContent) == c("OrganismName", "Group"))

library(ggfortify)
library(ggplot2)

WholeGenomeContentPCA <- prcomp(WholeGenomeContent[, -Species], scale = TRUE)
print(WholeGenomeContentPCA)
autoplot(WholeGenomeContentPCA, data=WholeGenomeContent, size=0.5, label=TRUE, label.size=1, loadings=TRUE, loadings.label=TRUE, label.vjust=1.5, colour='Group')

FourFoldContentPCA <- prcomp(FourFoldContent[, -SpeciesF], scale = TRUE)
print(FourFoldContentPCA)
autoplot(FourFoldContentPCA, data=FourFoldContent, size=0.5, label=TRUE, label.size=1, loadings=TRUE, loadings.label=TRUE, label.vjust=1.5, colour='grey53')

ff_pc <- data.frame(FourFoldContentPCA[['x']])
write.csv(ff_pc, '/home/emulciber/MitoClub/devil-worm/Body/3Results/11_PCA_FF_projections.csv', sep=';')


# а теперь с тандемными повторами
TRF <- read.table('/home/emulciber/MitoClub/devil-worm/Body/3Results/13_TRFoutputAnalysis.csv', head = TRUE, sep = ';')
TRF <- na.omit(TRF)
TRF[c(1,3:7)] <- list(NULL)
rownames(TRF) <- TRF$Species
SpeciesT <- which(colnames(TRF) == c("Species", "Group"))
TRF_PCA <- prcomp(TRF[, -SpeciesT], scale = TRUE)
print(TRF_PCA)
autoplot(TRF_PCA, data=TRF, size=0.5, label=TRUE, label.size=1, loadings=TRUE, loadings.label=TRUE, label.vjust=1.5, colour='Group')

trf_pc <- data.frame(TRF_PCA[["x"]])
write.csv(trf_pc, '/home/emulciber/MitoClub/devil-worm/Body/3Results/11_PCA_TRF_projections.csv', sep=';')

# и fourfold с тандемными повторами
FF_TRF <- FourFoldContent
names(FF_TRF)[names(FF_TRF) == 'OrganismName'] <- 'Species'
FF_TRF <- merge(x = FF_TRF, y = TRF, by = c("Species", 'Group'), all.x = TRUE)
rownames(FF_TRF) <- FF_TRF$Species
SpeciesFT <- which(colnames(FF_TRF) == c("Species", "Group"))
FF_TRF <- na.omit(FF_TRF)
FF_TRF[c(7:10)] <- list(NULL)

FF_TRF_PCA <- prcomp(FF_TRF[, -SpeciesFT], scale = TRUE)
print(FF_TRF_PCA)
autoplot(FF_TRF_PCA, data=FF_TRF, size=0.5, label=TRUE, label.size=1, loadings=TRUE, loadings.label=TRUE, label.vjust=1.5, colour='grey53')

# по 3 позиции
TrdPosContent <- read.table('/home/emulciber/MitoClub/devil-worm/Body/3Results/09_3rd_Pos_Content.csv', head = TRUE, sep = ';')

TrdPosContent[c(2:65)] <- list(NULL)

rownames(TrdPosContent) <- TrdPosContent$OrganismName

SpeciesT <- which(colnames(TrdPosContent) == c("OrganismName", "Group"))

TrdPosContentPCA <- prcomp(TrdPosContent[, -SpeciesT], scale = TRUE)
print(TrdPosContentPCA)
autoplot(TrdPosContentPCA, data=TrdPosContent, size=0.5, label=TRUE, label.size=1, loadings=TRUE, loadings.label=TRUE, label.vjust=1.5, colour='grey53')

tp_pc <- data.frame(TrdPosContentPCA[['x']])
write.csv(tp_pc, '/home/emulciber/MitoClub/devil-worm/Body/3Results/11_PCA_3rd_Pos_projections.csv', sep=';')
