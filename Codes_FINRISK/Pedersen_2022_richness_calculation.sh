#!/bin/bash
#SBATCH --job-name=richness
#SBATCH --account=project_2005073
#SBATCH --time=12:00:00
#SBATCH --partition=small
#SBATCH --ntasks=40
#SBATCH --mem-per-cpu=1500
##SBATCH --mail-type=BEGIN #uncomment to enable mail

get_seeded_random()
{
  seed="$1"
  openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt \
    </dev/zero 2>/dev/null
}

for name in $(cat "/scratch/project_2005073/USER_WORKSPACES/matti/Pedersen_2022/FR02_samplenames.txt")
do
while [`jobs | wc -l` -ge 40 ]
do
sleep 5
done
( nrow=$(wc -l "/scratch/project_2005073/USER_WORKSPACES/matti/Pedersen_2022/sorted/${name}.csv" | awk '{print $1}')
if [ $nrow -ge 100000 ]
then
echo "Processing sample $name"
for i in {1..5}
do
seed=$(( 42 * $i ))
sort --random-source=<(get_seeded_random ${seed}) -R /scratch/project_2005073/USER_WORKSPACES/matti/Pedersen_2022/sorted/${name}.csv | head -n 100000 | awk '{print $2}' | sort -u | wc -l | awk -v awkname="$name" '{print awkname, $0}' >> /scratch/project_2005073/USER_WORKSPACES/matti/Pedersen_2022/richness/${name}_richness.txt
done
fi ) &
done
wait

cat /scratch/project_2005073/USER_WORKSPACES/matti/Pedersen_2022/richness/*_richness.txt > /scratch/project_2005073/USER_WORKSPACES/matti/Pedersen_2022/richness_5.txt
awk '{seen[$1]+=$2; count[$1]++} END{for (x in seen)print x, seen[x]/count[x]}' /scratch/project_2005073/USER_WORKSPACES/matti/Pedersen_2022/richness_5.txt > /scratch/project_2005073/USER_WORKSPACES/matti/Pedersen_2022/richness_means.txt
