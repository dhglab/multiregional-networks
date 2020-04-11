# Shell script to run ARACNE-AP
# Gokul Ramaswami 10-10-2016
# DEPENDENCY: module load java/1.8.0_77

#$ -N ARACNE 
#$ -cwd 
#$ -l h_data=10G,h_rt=12:00:00,highp 
#$ -V 

# specify directories for data, output, software, and figures
DATDIR="../data/"
OUTDIR="../processed_data/"
FIGDIR="../figures/"
SOFTDIR="../software/"

TISSUE=$1 # Take tissue from command line
# Check that user specified tissue
if [ -z "$TISSUE" ]
then
	echo Please provide the tissue as input
	exit
fi

mkdir $OUTDIR$TISSUE
mkdir $FIGDIR$TISSUE
# calculate Threshold
java -Xmx5G -jar $SOFTDIR\ARACNE.jar -e $DATDIR\Expr_$TISSUE.txt -o $OUTDIR\$TISSUE --tfs $DATDIR\TFs_totalBrain.txt --seed 1 --pvalue 1E-8 --calculateThreshold

# Bootstrap networks
for i in {1..100}
do
java -Xmx5G -jar $SOFTDIR\ARACNE.jar -e $DATDIR\Expr_$TISSUE.txt -o $OUTDIR$TISSUE --tfs DATDIR\TFs_totalBrain.txt --pvalue 1E-8 --seed $i
done

# consensus from bootstraps
java -Xmx5G -jar $SOFTDIR\ARACNE.jar -o $OUTDIR$TISSUE --consolidate