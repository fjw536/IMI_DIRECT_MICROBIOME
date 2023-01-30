#!/bin/bash
#SBATCH --job-name=bwa-mem2_missing_3
#SBATCH --account=project_2005073
#SBATCH --time=24:00:00
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=100G
##SBATCH --mail-type=BEGIN #uncomment to enable mail

module load bwa-mem2/2.2
module load samtools/1.14

round=${SLURM_ARRAY_TASK_ID}

for i in {1..27}
	do
	# set the input file to process
	name=$(sed "${i}q;d" "/scratch/project_2005073/USER_WORKSPACES/matti/Pedersen_2022/FR02_samplenames.txt")
	#map with -r 1 and -D 0.3 to the gene database
	bwa-mem2 mem -t 8 -r 1 -D 0.3 /scratch/project_2005073/DATABASES/MDB/IGC.fa /scratch/project_2005073/DATA/${name}.R1.trimmed.filtered.fastq.gz /scratch/project_2005073/DATA/${name}.R2.trimmed.filtered.fastq.gz > /scratch/project_2005073/USER_WORKSPACES/matti/Pedersen_2022/${name}_out.sam

	#fix mate reads and sort
	samtools fixmate -m /scratch/project_2005073/USER_WORKSPACES/matti/Pedersen_2022/${name}_out.sam /scratch/project_2005073/USER_WORKSPACES/matti/Pedersen_2022/${name}_out_fixmate.sam
	samtools sort -o /scratch/project_2005073/USER_WORKSPACES/matti/Pedersen_2022/${name}_out_sorted.sam /scratch/project_2005073/USER_WORKSPACES/matti/Pedersen_2022/${name}_out_fixmate.sam
	#remove duplicates
	samtools markdup -r /scratch/project_2005073/USER_WORKSPACES/matti/Pedersen_2022/${name}_out_sorted.sam /scratch/project_2005073/USER_WORKSPACES/matti/Pedersen_2022/${name}_out_nodup.sam
	#filter to only mapped reads, with >= 95% identity (according to edit distance) across a window of >= 75 bases, and with MAPQ >= 20
	#a hit is counted for a read pair if both of the reads in a pair map to the same gene, and at least one of these reads surpasses the specified criteria
	#however, multimapped read pairs are completely removed
	#the results per sample are saved in a table with combinations of reads and genes on each row
	samtools view -F 0x4 -q 20 /scratch/project_2005073/USER_WORKSPACES/matti/Pedersen_2022/${name}_out_nodup.sam | awk '{match($12, /NM:i:([0-9]*)/, nm)} $7 =="=" && length($10) >= 75 && nm[1] <= length($10)*0.05 {print $1, $3}' | sort -u | awk '{k=($1 FS $2)} {a[$1]++;b[$1]=k}END{for(x in a)if(a[x]==1)print b[x]}' > /scratch/project_2005073/USER_WORKSPACES/matti/Pedersen_2022/sorted/${name}.csv
	#intermediary files are removed
	rm /scratch/project_2005073/USER_WORKSPACES/matti/Pedersen_2022/${name}_out*
done
