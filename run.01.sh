# Get and map files

# Get chromosome names
  grep '>' /data/tunglab/shared/genomes/panubis1/Panubis1.0.fa > 00_chroms.txt; sed -i 's/>//g' 00_chroms.txt

# Pull files from SRA 
	sbatch --array=1-44 --mem=16G run01_getSRA.sh

#Trim & map 
	sbatch --array=1-44 --mem=28G run02_trimmap.sh 
