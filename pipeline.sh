#!/bin/bash
#Store fastq data from the server to local
mkdir fastq | cp -r /localdisk/data/BPSM/ICA1/fastq/* ./fastq

#perform fastqc, creat /fastq_output and output to this floder
mkdir fastq_output | fastqc ./fastq/*.fq.gz -o ./fastq_output

#generate the number and quality of sequence from fastqc,save a txt which including filename,totalSequences,poorQuality sequence,pass/fail

unzip ./fastq_output/"*.zip" -d ./fastq_output/

for subdir in $(find fastq_output -mindepth 1 -type d)
do

    if [ -f "$subdir/fastqc_data.txt" ]
    then
                filename=$(grep -P '^\Filename' "$subdir/fastqc_data.txt" | awk -F'\t' '{print $2}')
                totalSeq=$(grep -P '^\Total Sequences' "$subdir/fastqc_data.txt" | awk -F'\t' '{print $2}')
                poorQuality=$(grep -P '^\Sequences flagged as poor quality' "$subdir/fastqc_data.txt" | awk -F'\t' '{print $2}')
                status=$(grep -P '^>>Basic Statistics' "$subdir/fastqc_data.txt" | awk -F'\t' '{print $2}')

                echo -e "$filename\t$totalSeq\t$poorQuality\t$status">> ./fastqc_output.txt

    fi
done

#Store Tcongo_genome data from the server to loacl
mkdir Tcongo_genome | cp -r /localdisk/data/BPSM/ICA1/Tcongo_genome/* ./Tcongo_genome

#gunzip the file
gunzip ./Tcongo_genome/*

#bowtie-build
cd Tcongo_genome
bowtie2-build ./TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta genome

filepath="../fastq/Tco2.fqfiles"

mkdir -p genome_bowtie_output
mkdir -p genome_samtool_output
mkdir -p genome_samtool_sort
mkdir -p genome_bedtools_output
while read -r line; do
    samplename=$(echo "$line" | awk '{print $1}')
    sampletype=$(echo "$line" | awk '{print $2}')
    time=$(echo "$line" | awk '{print $4}')
    Treatment=$(echo "$line" | awk '{print $(5)}')
    end1=$(echo "$line" | awk '{print $(NF-1)}')
    end2=$(echo "$line" | awk '{print $(NF)}')

    echo "$samplename $sampletype $time $end1 $end2"
#create the file floder,and help calssification
    mkdir -p ./genome_bedtools_output/${sampletype}_${time}_${Treatment}
#bowtie
    bowtie2 -x genome -1 ../fastq/$end1 -2 ../fastq/$end2 -S ./genome_bowtie_output/${samplename}.sam
#samtools
    samtools view -bS ./genome_bowtie_output/${samplename}.sam > ./genome_samtool_output/${samplename}.bam
    samtools sort ./genome_samtool_output/${samplename}.bam -o ./genome_samtool_sort/${samplename}_sort.bam
#bedtools,and save the txt to the classified folder    
    bedtools coverage -counts -a /localdisk/data/BPSM/ICA1/TriTrypDB-46_TcongolenseIL3000_2019.bed -b ./genome_samtool_sort/${samplename}_sort.bam > ./genome_bedtools_output/${sampletype}_${time}_${Treatment}/${samplename}.txt

done < <(tail -n +2 "$filepath")

# first, paste together the genenames descriptions and gene counts of groups
parent_folder="./genome_bedtools_output/"
mkdir -p ../average_output

for subdir in "$parent_folder"/*
do
        i=1
        temp_file=$(mktemp)
        for file in "$subdir"/*
        do
            if [ "$i" -eq 1 ]
            then
                awk 'BEGIN{FS="\t"; OFS="\t";}{print $4,$5,$6;}' $file > $temp_file
                ((i++))
            else
                next_temp_file=$(mktemp)
                paste $temp_file <(awk -F'\t' '{print $6}' $file) > $next_temp_file
                mv $next_temp_file $temp_file
                ((i++))
            fi
        done
        cp $temp_file ../average_output/${subdir##*/}.txt
        rm $temp_file
#calulate the mean,and delete the processing file
        awk 'BEGIN{FS=OFS="\t"} {sum=0; for(i=3; i<=NF; i++) sum+=$i; avg=sum/(NF-2); print $1, $2, avg}' ../average_output/${subdir##*/}.txt > ../average_output/${subdir##*/}_average.txt
        rm ../average_output/${subdir##*/}.txt
done


mkdir -p ../final_output
cd ../average_output

# generate "fold change" data for the "group-wise" comparisons

#over time between clones
mkdir -p ../final_output/clone1
mkdir -p ../final_output/clone2

reference1="./WT_0_Uninduced_average.txt"
reference2="./WT_24_Induced_average.txt"
reference3="./WT_24_Uninduced_average.txt"
reference4="./WT_48_Induced_average.txt"
reference5="./WT_48_Uninduced_average.txt"

#because in reference groups,gene expression can be zero, and cannot divide,so use log2(A+1)-log2(b+1) to present log2(fold change)
paste "./Clone1_0_Uninduced_average.txt" "$reference1" | awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, (log($3+1)/log(2) - log($6+1)/log(2))}' > ../final_output/clone1/Clone1_0_Uninduced_fold_change.txt
paste "./Clone2_0_Uninduced_average.txt" "$reference1" | awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, (log($3+1)/log(2) - log($6+1)/log(2))}' > ../final_output/clone2/Clone2_0_Uninduced_fold_change.txt

paste "./Clone1_24_Induced_average.txt" "$reference2" | awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, (log($3+1)/log(2) - log($6+1)/log(2))}' > ../final_output/clone1/Clone1_24_Induced_fold_change.txt
paste "./Clone2_24_Induced_average.txt" "$reference2" | awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, (log($3+1)/log(2) - log($6+1)/log(2))}' > ../final_output/clone2/Clone2_24_Induced_fold_change.txt
paste "./Clone1_24_Uninduced_average.txt" "$reference3" | awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, (log($3+1)/log(2) - log($6+1)/log(2))}' > ../final_output/clone1/Clone1_24_Uninduced_fold_change.txt
paste "./Clone2_24_Uninduced_average.txt" "$reference3" | awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, (log($3+1)/log(2) - log($6+1)/log(2))}' > ../final_output/clone2/Clone2_24_Uninduced_fold_change.txt

paste "./Clone1_48_Induced_average.txt" "$reference4" | awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, (log($3+1)/log(2) - log($6+1)/log(2))}' > ../final_output/clone1/Clone1_48_Induced_fold_change.txt
paste "./Clone2_48_Induced_average.txt" "$reference4" | awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, (log($3+1)/log(2) - log($6+1)/log(2))}' > ../final_output/clone2/Clone2_48_Induced_fold_change.txt

paste "./Clone1_48_Uninduced_average.txt" "$reference5" | awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, (log($3+1)/log(2) - log($6+1)/log(2))}' > ../final_output/clone1/Clone1_48_Uninduced_fold_change.txt
paste "./Clone2_48_Uninduced_average.txt" "$reference5" | awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, (log($3+1)/log(2) - log($6+1)/log(2))}' > ../final_output/clone2/Clone2_48_Uninduced_fold_change.txt


#sort
folder_path="../final_output/clone2"
for file in "$folder_path"/*
do
        sort -k4,4nr -t$'\t' <(awk 'function abs(x){return x<0?-x:x} BEGIN {FS="\t"; OFS="\t"} {print $1, $2, $3, abs($NF)}' $file | sort -k4,4nr)| cut -f1-3 > "${file%.*}.sorted.txt"
        rm "${file}"
done

folder_path="../final_output/clone1"
for file in "$folder_path"/*
do
        sort -k4,4nr -t$'\t' <(awk 'function abs(x){return x<0?-x:x} BEGIN {FS="\t"; OFS="\t"} {print $1, $2, $3, abs($NF)}' $file | sort -k4,4nr)| cut -f1-3 > "${file%.*}.sorted.txt"
        rm "${file}"
done

echo "Well done! That's all!"
