# Get genomic contexts that matter

# genes, enhancers, and cpg islands from Vilgalys & Fogel et al. bioRxiv are used to define regions of the genome. 
## for genes we only want protein coding genes
mkdir contexts 
cp ../github/VilgalysFogel_Amboseli_admixture/03_baboon-genomic-resources/Resources/enhancers_lifted_to_baboon.bed ./contexts/
cp ../github/VilgalysFogel_Amboseli_admixture/03_baboon-genomic-resources/Resources/cpg_islands.panubis1.emboss.bed ./contexts/
module load Anaconda3/2019.10-gcb02; conda activate r4-base; R
library(data.table); all <- fread("../my_genomes/panubis1/Panubis1_gtf_readable.txt", fill=T)
exon <- all[all$V3 == 'exon',]
gene <- all[all$V3 == 'gene',]; gene <- gene[gene$V13 == 'protein_coding"',]
exon <- exon[exon$V13 %in% gene$V12,]
gene$V1 <- paste("chr", gene$V1, sep="")
exon$V1 <- paste("chr", exon$V1, sep="")
exon <- exon[exon$V12 == 'mRNA',]
write.table(gene[,c(1,4:5,9,10,13)], "./contexts/genes.bed", row.names=F, col.names=F, sep="\t", quote=F)
write.table(exon[,c(1,4:5,9,11,10,16,12)], "./contexts/exons.bed", row.names=F, col.names=F, sep="\t", quote=F)
g_for <- gene[gene$V7 == "+",]
g_rev <- gene[gene$V7 == "-",]
g_for$start <- g_for$V4-2000; g_for$start[g_for$start <= 0] <- 1; g_for$end <- g_for$V4 
g_rev$start <- g_rev$V5; g_rev$end <- g_rev$V5+2000
rbind(g_for, g_rev) -> prom
write.table(prom[,c(1,17:18, 9,10,13)], "./contexts/prom_2kb.bed", row.names=F, col.names=F, sep="\t", quote=F)
cpg <- fread("./contexts/cpg_islands.panubis1.emboss.bed", header=F) 
tmp1 <- tmp2 <- cpg[cpg$V1 %like% 'chr',]
tmp1$V3 <- tmp1$V2-1; tmp1$V2 <- tmp1$V3 - 2000; tmp1$V2[tmp1$V2 <= 0] <- 1
tmp2$V2 <- tmp2$V3+1; tmp2$V3 <- tmp2$V2 + 2000
shore <- rbind(tmp1, tmp2)
write.table(shore, "./contexts/cpg_shore_2kb.bed", row.names=F, col.names=F, sep="\t", quote=F)
load("n44.raw_count_data.RData")
write.table(info[,1:4], "./sites.bed", row.names=F, col.names=F, sep="\t", quote=F)
quit(save='no')
conda deactivate 
module load bedtools2 

bedtools intersect -a sites.bed -b ./contexts/exons.bed -u > sites.exon.bed
bedtools intersect -a sites.bed -b ./contexts/genes.bed -u > sites.gene.bed
bedtools intersect -a sites.bed -b ./contexts/prom_2kb.bed -u > sites.promoter.bed
bedtools intersect -a sites.bed -b ./contexts/cpg_islands.panubis1.emboss.bed -u > sites.island.bed
bedtools intersect -a sites.bed -b ./contexts/cpg_shore_2kb.bed -u > sites.shore.bed
bedtools intersect -a sites.bed -b ./contexts/enhancers_lifted_to_baboon.bed -u > sites.enhancer.bed
bedtools intersect -a sites.bed -b ./contexts/*.bed -v > sites.unannotated.bed

bedtools intersect -a ./contexts/exons.bed  -b sites.bed -u > featured.exon.bed
bedtools intersect -a ./contexts/genes.bed  -b sites.bed -u > featured.gene.bed
bedtools intersect -a ./contexts/prom_2kb.bed  -b sites.bed -u > featured.promoter.bed
bedtools intersect -a ./contexts/cpg_islands.panubis1.emboss.bed  -b sites.bed -u > featured.island.bed
bedtools intersect -a ./contexts/cpg_shore_2kb.bed  -b sites.bed -u > featured.shore.bed
bedtools intersect -a ./contexts/enhancers_lifted_to_baboon.bed  -b sites.bed -u > featured.enhancer.bed

# Integrate into our data frame 
conda activate r4-base; R
library(data.table); load("n44.raw_count_data.RData")
info$in.exon <- info$in.gene <- info$in.gene_not_exon <- info$in.promoter <- info$in.island <- info$in.shore <- info$in.enhancer <- info$in.unannotated <- 0
fread("sites.exon.bed") -> tmp; info$in.exon[info$site %in% tmp$V4] <- 1
fread("sites.gene.bed") -> tmp; info$in.gene[info$site %in% tmp$V4] <- 1
info$in.gene_not_exon[info$in.gene == 1 & info$in.exon == 0] <- 1
fread("sites.promoter.bed") -> tmp; info$in.promoter[info$site %in% tmp$V4] <- 1
fread("sites.island.bed") -> tmp; info$in.island[info$site %in% tmp$V4] <- 1
fread("sites.shore.bed") -> tmp; info$in.shore[info$site %in% tmp$V4] <- 1
fread("sites.enhancer.bed") -> tmp; info$in.enhancer[info$site %in% tmp$V4] <- 1
fread("sites.unannotated.bed") -> tmp; info$in.unannotated[info$site %in% tmp$V4] <- 1
save.image("n44.raw_count_data.RData")


# regions fitting the concensus baboon phylogeny will be identified using baboons from the Baboon Genome Diversity Panel (Rogers et al. 2019), remapped to the panubis1 genome in Vilgalys & Fogel et al. 

# regions are identified in ./CM/alternate_phylogeny/

# get vcf files 
sbatch --array=1-20 --mem=16G run.01.get_unadmixed.sh
# returns filtered vcf file per chromosome, using genotypes already called and the list of files samples in 00_allspecies.list 

module load bcftools; bcftools concat -O z -o BGDP_samples.vcf.gz ./chrom_vcfs/02.*gz
module load tabix; tabix BGDP_samples.vcf.gz

module load Anaconda3/2019.10-gcb02; conda activate r4-base; R


# get partial vcf files. 
## 1 Mb windows, sliding every 100kb 
for f in `seq 1 100000 218172882`; do tabix BGDP_samples.vcf.gz 1:$f-$(($f+1000000)) > partial_vcfs/chr1.$f.vcf; done
for f in `seq 1 100000 193660750`; do tabix BGDP_samples.vcf.gz 2:$f-$(($f+1000000)) > partial_vcfs/chr2.$f.vcf; done
for f in `seq 1 100000 184919515`; do tabix BGDP_samples.vcf.gz 3:$f-$(($f+1000000)) > partial_vcfs/chr3.$f.vcf; done
for f in `seq 1 100000 182120902`; do tabix BGDP_samples.vcf.gz 4:$f-$(($f+1000000)) > partial_vcfs/chr4.$f.vcf; done
for f in `seq 1 100000 173900761`; do tabix BGDP_samples.vcf.gz 5:$f-$(($f+1000000)) > partial_vcfs/chr5.$f.vcf; done
for f in `seq 1 100000 167138247`; do tabix BGDP_samples.vcf.gz 6:$f-$(($f+1000000)) > partial_vcfs/chr6.$f.vcf; done
for f in `seq 1 100000 161768468`; do tabix BGDP_samples.vcf.gz 7:$f-$(($f+1000000)) > partial_vcfs/chr7.$f.vcf; done
for f in `seq 1 100000 140274886`; do tabix BGDP_samples.vcf.gz 8:$f-$(($f+1000000)) > partial_vcfs/chr8.$f.vcf; done
for f in `seq 1 100000 127591819`; do tabix BGDP_samples.vcf.gz 9:$f-$(($f+1000000)) > partial_vcfs/chr9.$f.vcf; done
for f in `seq 1 100000 126462689`; do tabix BGDP_samples.vcf.gz 10:$f-$(($f+1000000)) > partial_vcfs/chr10.$f.vcf; done
for f in `seq 1 100000 125913696`; do tabix BGDP_samples.vcf.gz 11:$f-$(($f+1000000)) > partial_vcfs/chr11.$f.vcf; done
for f in `seq 1 100000 123343450`; do tabix BGDP_samples.vcf.gz 12:$f-$(($f+1000000)) > partial_vcfs/chr12.$f.vcf; done
for f in `seq 1 100000 106849001`; do tabix BGDP_samples.vcf.gz 13:$f-$(($f+1000000)) > partial_vcfs/chr13.$f.vcf; done
for f in `seq 1 100000 106654974`; do tabix BGDP_samples.vcf.gz 14:$f-$(($f+1000000)) > partial_vcfs/chr14.$f.vcf; done
for f in `seq 1 100000 91985775`; do tabix BGDP_samples.vcf.gz 15:$f-$(($f+1000000)) > partial_vcfs/chr15.$f.vcf; done
for f in `seq 1 100000 91184193`; do tabix BGDP_samples.vcf.gz 16:$f-$(($f+1000000)) > partial_vcfs/chr16.$f.vcf; done
for f in `seq 1 100000 74525926`; do tabix BGDP_samples.vcf.gz 17:$f-$(($f+1000000)) > partial_vcfs/chr17.$f.vcf; done
for f in `seq 1 100000 72894408`; do tabix BGDP_samples.vcf.gz 18:$f-$(($f+1000000)) > partial_vcfs/chr18.$f.vcf; done
for f in `seq 1 100000 72123344`; do tabix BGDP_samples.vcf.gz 19:$f-$(($f+1000000)) > partial_vcfs/chr19.$f.vcf; done
for f in `seq 1 100000 50021108`; do tabix BGDP_samples.vcf.gz 20:$f-$(($f+1000000)) > partial_vcfs/chr20.$f.vcf; done


## in R 
library(SNPRelate); library(ctc); library(ape)

# make GDS object and get phylogenetic tree for all samples 
snpgdsVCF2GDS("BGDP_samples.vcf.gz", "full.gds")
snpgdsSummary("full.gds")
genofile <- snpgdsOpen("full.gds")
dissMatrix  <-  snpgdsDiss(genofile , sample.id=NULL, snp.id=NULL, autosome.only=TRUE,remove.monosnp=TRUE, maf=NaN, missing.rate=NaN, num.thread=10, verbose=TRUE)
hclust(as.dist(dissMatrix$diss)) -> tree
tree$labels <- c('Anubis', 'Anubis', 'Anubis', 'Anubis', 'Yellow', 'Yellow', 'Hamadryas', 'Hamadryas', 'Kinda', 'Kinda', 'Kinda', 'Guinea', 'Guinea', 'Chacma', 'Chacma')
# note the first panubis individual (30877) is the one that seems to contain a lot of yellow ancestry, let's not use him for anything

#hc2Newick(tree) -> Newick
vcv.phylo(as.phylo(tree))

