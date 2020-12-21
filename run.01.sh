# Get and map files

# Get chromosome names
  grep '>' /data/tunglab/shared/genomes/panubis1/Panubis1.0.fa > 00_chroms.txt; sed -i 's/>//g' 00_chroms.txt

# Pull files from SRA 
  sbatch --array=1-44 --mem=16G run01_getSRA.sh

#Trim & map (this needs ~25gigs of memory per sample)
	# NOTE 1: you need to call the methratio.py script here, this can be downloaded from the BSMAP website. Also, methratio.py only works on SAM files, so the mapping output is set to SAM rather than BAM. 
	# NOTE 2: map to the baboon and the lambda phage genomes separately. Not doing here as I'm just keeping the previously reported conversion rates. 
  sbatch --array=1-44 --mem=28G run02_trimmap.sh 

# Quality control
	# Not done for panubis1
	#Bisulfite conversion rate 
  for f in `ls *methratio.lam.txt`; do $f; grep 'lambda' $f | awk '{ SUM += $8} END { print SUM}'; grep 'lambda' $f | awk '{ SUM += $7} END { print SUM}'; done
	#Msp1 Digest
  for k in `ls *panu.sam`; do $k; grep -P '\t(CGG|TGG)' $k | wc -l ; 
	grep -P '(CCA|CCG)\t' $k | wc -l ; 
	grep -P '(CGG|TGG)\t' $k | wc -l ; 
	grep -P '\t(CCA|CCG)' $k | wc -l ;
  done
