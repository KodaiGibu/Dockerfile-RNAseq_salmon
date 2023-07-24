#!/bin/sh

#salmonで使用するindexのパスを事前に代入する
db="/mnt/Adig/Adig"

#Rawデータのstat
# リード数
for file in `ls|grep gz`
do
seqkit -j 2 stat $file >> Sum1.txt
done
grep -v "file" Sum1.txt | awk '{print $1 $4}' > Sum1r.txt

#trimmomaticを用いたアダプタートリミング
#サンプルリストの作成
ls *_1.fq.gz | sed -e 's/_1.fq.gz//g' > list.txt

#トリミング＆フィルタリング
for file in `cat list.txt`; do
TrimmomaticPE -threads 2 -phred33 -trimlog log.txt ${file}_1.fq.gz ${file}_2.fq.gz ${file}_R1_001_trimmed.fastq.gz ${file}_R1_001_unpaired.fastq.gz ${file}_R2_001_trimmed.fastq.gz ${file}_R2_001_unpaired.fastq.gz ILLUMINACLIP:adapters.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:50
done

#マッピングしなかったファイルの除去。このファイルは出力しなくてもいい。
rm *unp*

#トリミングしたファイルを一つのフォルダにまとめる。
mkdir trimmed
mv *trimmed.fastq.gz ./trimmed

#fastqcで処理したファイルを格納するフォルダを作成
mkdir fastqc
for file in `cat list.txt`; do
fastqc -t 10 --nogroup -o fastqc -f fastq ./trimmed/${file}_R1_001_trimmed.fastq.gz ./trimmed/${file}_R2_001_trimmed.fastq.gz
done

#indexの作成
#singularity exec /usr/local/biotools/s/salmon\:1.9.0--h7e5ed60_1 salmon index -t ./rna.fna -i ./Adig -k 31

mkdir salmon_result

#定量化
for seqlib in `cat list.txt`; do
result_dir=${seqlib}_exp_salmon
salmon quant -i ${db} -p 2 --validateMappings -l A -1 ./trimmed/${seqlib}_R1_001_trimmed.fastq.gz -2 ./trimmed/${seqlib}_R2_001_trimmed.fastq.gz -o ./salmon_result/${result_dir}
done

#マッピング率を抽出し、一つのファイルにまとめる
for file in `cat list.txt`; do
grep "Mapping rate" salmon_result/${file}_exp_salmon/logs/salmon_quant.log | sed 's/.*= //' >> MappingRate.txt
done

#list.txtと結合してラベルを付ける
paste list.txt MappingRate.txt > MappingRate2.txt

#中間ファイルの削除。リネーム。
rm MappingRate.txt
mv MappingRate2.txt MappingRate.txt
