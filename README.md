# RNAseq_salmon
Dockerfile for RNA-seq with salmon

# 実践編
実践編では、オニヒトデのRNA-seqデータ（鎌田さんサンプル）を用いて行う。
解析方法は<a href="https://qiita.com/gibukod/private/5540a8fc1eaab88b6bba">salmonを用いる方法</a>で実行する。

<a href="https://www.dropbox.com/scl/fo/j506ohieduz6v7pqv09al/h?rlkey=33tpms25x00xkgapyuot10mf0&dl=0">raw dataおよびスクリプトはここからダウンロード</a>
遺伝研スパコンで、すでにインストールされたソフトを使う場合はcode02.shを使うこと。



## 必要なソフトのインストール
condaによるインストールがシンプルかつ簡単。
定量化ソフトであるrsemはcondaでインストールできるバージョンがかなり古いので、githubからソースコードをダウンロードし、コンパイルした。
(condaでのrsemのインストールは可能だが、バージョンが古く、今回用いるコードでは動かない。古いバージョンのコードでの書き直しが必要。）
```sh:ソフトのダウンロード
#condaによるインストール。
#ソフトのインストール。
conda create -n rnaseq_mapping
conda activate rnaseq_mapping

conda install -c bioconda trimmomatic -y
conda install -c bioconda fastqc -y
conda install -c conda-forge seqkit -y

conda update -c conda-forge trimmomatic -y
conda update -c conda-forge fastqc -y
conda update -c conda-forge seqkit -y

#salmonはcondaでインストールするとエラーがでた（バージョンも古い）。
#バイナリファイルをダウンロードしてpathを通す。
wget https://github.com/COMBINE-lab/salmon/releases/download/v1.10.0/salmon-1.10.0_linux_x86_64.tar.gz
tar -xvf salmon-1.10.0_linux_x86_64.tar.gz

#salmonの最新版はv1.10.1だが（2023/05/25時点）コンパイルの必要があるので、今回はv1.10.0をダウンロードした。

```
## リファレンスRNAデータのダウンロードおよびindexの作成
リファレンスとなるトランスクリプトームデータをダウンロードする。
今回は以下のサイトよりダウンロード。
<p><a href="https://bioinfo.szn.it/acanthaster-planci/">BAC: HPC AND BIOINFORMATICS AT SZN</a></p>

AplaはAcanthaster planciの略。

```sh:リファレンスの準備
#リファレンスを入れるディレクトリを作成
mkdir db
cd db

#リファレンストランスクリプトームのfastaをダウンロード
wget http://bioinfo.szn.it/genoma-repo/Aplanci/REFSEQ/GCF_001949145.1_OKI-Apl_1.0_rna.fna.tbz2

#ダウンロードしたfastaを解凍
tar -xvf GCF_001949145.1_OKI-Apl_1.0_rna.fna.tbz2

#解凍前のファイルを削除
rm GCF_001949145.1_OKI-Apl_1.0_rna.fna.tbz2

#salmonでリファレンスとして利用するため、indexを作成
salmon index -t ./GCF_001949145.1_OKI-Apl_1.0_rna.fna -i ./Apla -k 31
#遺伝研スパコンの場合は以下
singularity exec /usr/local/biotools/s/salmon\:1.9.0--h7e5ed60_1 salmon index -t ./GCF_001949145.1_OKI-Apl_1.0_rna.fna -i ./Apla -k 31

#作成されたindexディレクトリを確認
ls Apla/
```
問題なくindexの作成が終了すると、Aplaというディレクトリが作成され、中にcomplete_ref_lens.binなどのファイルが生成されている。

## salmonによる解析の半自動化スクリプト
trimmomaticによるトリミング、fastqcによるクオリティチェック、salmonによるリードカウントを行うスクリプトを以下からダウンロード可能。
ペアエンド用。シングルエンドのfastqを使う場合はtrimmomaticのオプションを書き換える必要がある。

<a href="https://www.dropbox.com/scl/fi/0mr2p77uq1twryxla34tb/code01.sh?rlkey=znkjev1ivb8n3m8ylcfjp7hn1&dl=0">半自動化スクリプト</a>
<a href="https://www.dropbox.com/scl/fi/wla73v7kbm8jca0trjygs/code02.sh?rlkey=1x41a7wdh1fm1uesotk4laald&dl=0">半自動化スクリプト_遺伝研スパコンver</a>

<br>

ディレクトリにfastqファイル、adapters.fa、半自動化スクリプト、salmon用index（上記で作成）を置く。
```ディレクトリ内のファイルの確認
[gibukodai@at139 Ap]$ ls
St19_1.fq.gz
St19_2.fq.gz
St20_1.fq.gz
St20_2.fq.gz
St21_1.fq.gz
St21_2.fq.gz
St22_1_1.fq.gz
St22_1_2.fq.gz
adapters.fa
code01.sh
db
```
adapters.faはアダプター配列をfataファイルとして保存たものである。
```adapters.fa
>Read_1
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
>Read_2
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
```
<br>

スクリプト内の`db="/home/gibukodai/216/Ap/db/Apla"`を自分のsalmon用indexのAplaのパスに書き換える。

<br>
以下のようにスクリプトを実行し、解析が終わるまで待つ。

```sh:スクリプトの実行
#以下のコードによりバックグラウンドでジョブを実行
nohup bash code01.sh &
#実行を確認
top

#遺伝研スパコンの場合はqloginしてqsubで実行
qlogin
qsub code02.sh
#実行を確認。statがqwで待機中、rで実行中、Eqwでエラー。
qstat
#Eqwが出た場合は以下でジョブをデリート
qdel job_ID


#スクリプトに許可がない場合は実行できないので、以下のようにして許可を与える(permission deniedと表示される)
chmod +x code01.sh

```

<br>

## スクリプトの内容
### indexのpathを指定
ここで最初に作成したindexのpathを指定している。
解析前にスクリプトを編集して自身のindexのpathを書き込む必要がある。
```sh:indexのpath
#index path
db="/home/gibukodai/216/Ap/db/Apla"
```

<br>

### seqkitで初期ステータスの確認
seqkitでリード数などのデータの初期ステータスを抽出する。

```sh:初期値の記録
#Rawデータのstat
for file in `ls|grep gz`
do
seqkit -j 8 stat $file >> Sum1.txt
done
grep -v "file" Sum1.txt | awk '{print $1 $4}' > Sum1r.txt

#遺伝研スパコンの場合は以下
for file in `ls|grep gz`
do
singularity exec /usr/local/biotools/s/seqkit\:2.0.0--h9ee0642_0 seqkit -j 2 stat $file >> Sum1.txt
done
grep -v "file" Sum1.txt | awk '{print $1 $4}' > Sum1r.txt
```
<br>

### trimomaticによるアダプタートリミングおよび低クオリティリードの除去
trimmomaticを用いてアダプタートリミングおよび低クオリティリードの除去を一括で行う。この過程はfastpやcutadaptなど別のソフトを使用してもよい。

:::note warn
注意
サンプルリスト作成時のgzファイルの末尾はシーケンサによって異なるので、自分の使用するサンプルに合わせて書き換える必要がある（_1.fq.gzの部分）。
trimmomaticを使う場合は事前にadapters.faを作成し、アダプター配列を指定する必要がある。
:::

以下のコードでは、ペアエンドリードのサンプル名のみをリストに書き出し、そのリストを参照してfor文を用いて一気に処理する。

```sh:trimmomatic
#trimmomaticを用いたアダプタートリミング
#サンプルリストの作成
ls *_1.fq.gz | sed -e 's/_1.fq.gz//g' > list.txt

#trimmomaticの実行
for file in `cat list.txt`; do
trimmomatic PE -threads 4 -phred33 -trimlog log.txt\ 
${file}_1.fq.gz ${file}_2.fq.gz\ 
${file}_R1_001_trimmed.fastq.gz ${file}_R1_001_unpaired.fastq.gz\ 
${file}_R2_001_trimmed.fastq.gz ${file}_R2_001_unpaired.fastq.gz \
ILLUMINACLIP:adapters.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:50
done
#遺伝研スパコンの場合は以下
for file in `cat list.txt`; do
singularity exec /usr/local/biotools/t/trimmomatic\:0.39--hdfd78af_2 trimmomatic PE -threads 2 -phred33 -trimlog log.txt ${file}_1.fq.gz ${file}_2.fq.gz ${file}_R1_001_trimmed.fastq.gz ${file}_R1_001_unpaired.fastq.gz ${file}_R2_001_trimmed.fastq.gz ${file}_R2_001_unpaired.fastq.gz ILLUMINACLIP:adapters.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:50
done

#マッピングしなかったファイルの除去。このファイルは出力しなくてもいい。
rm *unp*

#トリミングしたファイルを一つのフォルダにまとめる。
mkdir trimmed
mv *trimmed.fastq.gz ./trimmed

```

<br>


### クリーニングしたリードのクオリティチェック
fastqcを用いてリードの処理し、出力されるhtmlファイルを開くことでクオリティを確認できる。
著者はいちいちスパコンからhtmlファイルをダウンロードしているが、vscodeなどのプラグインを利用することでssh接続したPC内でhtmlファイルを開くことができる。
トリミングで用いたリストを利用し、for文で一気に処理する。
```sh:クオリティチェック
#fastqcで処理したファイルを格納するフォルダを作成
mkdir fastqc

#fastqcの実行
for file in `cat list.txt`; do
fastqc -t 10 --nogroup -o fastqc -f fastq ./trimmed/${file}_R1_001_trimmed.fastq.gz ./trimmed/${file}_R2_001_trimmed.fastq.gz
done

#遺伝研スパコンの場合は以下
mkdir fastqc
for file in `cat list.txt`; do
singularity exec /usr/local/biotools/f/fastqc\:0.12.1--hdfd78af_0 fastqc -t 10 --nogroup -o fastqc -f fastq ./trimmed/${file}_R1_001_trimmed.fastq.gz ./trimmed/${file}_R2_001_trimmed.fastq.gz
done
```

<br>

### salmonを用いた転写産物の定量
salmonはリファレンス配列にリードを疑似アライメントをすることで定量化を行っている。

```sh:転写産物の定量
#salmonの実行結果を格納するディレクトリを作成
mkdir salmon_result

#salmonの実行
for seqlib in `cat list.txt`; do
result_dir=${seqlib}_exp_salmon
salmon quant -i ../ref/transcripts_index_salmon -p 4 --validateMappings -l A -1 ${seqlib}_R1_001_trimmed.fastq.gz -2 ${seqlib}_R2_001_trimmed.fastq.gz -o ./${result_dir}
done

#遺伝研スパコンの場合は以下
for seqlib in `cat list.txt`; do
result_dir=${seqlib}_exp_salmon
singularity exec /usr/local/biotools/s/salmon\:1.9.0--h7e5ed60_1 salmon quant -i ${db} -p 2 --validateMappings -l A -1 ./trimmed/${seqlib}_R1_001_trimmed.fastq.gz -2 ./trimmed/${seqlib}_R2_001_trimmed.fastq.gz -o ./salmon_result/${result_dir}
done

```
うまく実行されれば、サンプル名_exp_salmonというフォルダにカウントデータが生成され、生成されたquant.sfにというファイルにカウントデータが出力されている。
<br>
以下のコードを実行することでlogからマッピング率を抽出している。
```sh:マッピング率の抽出
#マッピング率を抽出し、一つのファイルにまとめる
for file in `cat list.txt`; do
grep "Mapping rate" salmon_result/${file}_exp_salmon/logs/salmon_quant.log | sed 's/.*= //' >> MappingRate.txt
done

#list.txtと結合してラベルを付ける
paste list.txt MappingRate.txt > MappingRate2.txt

#中間ファイルの削除。リネーム。
rm MappingRate.txt
mv MappingRate2.txt MappingRate.txt

```
<br>


### 出力ファイル
実行が終わると以下のようにファイルが出力されている。
trimmedにトリミング済みfastq、fastqcにクオリティの記されたhtmlファイル、salmon_resultにリードカウントの結果が格納されている。
次の解析にはsalmon_resultを用いる。

```:マッピング率の抽出
[gibukodai@at138 Ap]$ ls
St19_1.fq.gz
St19_2.fq.gz
St20_1.fq.gz
St20_2.fq.gz
St21_1.fq.gz
St21_2.fq.gz
St22_1_1.fq.gz
St22_1_2.fq.gz
Sum1.txt 
Sum1r.txt
code02.sh
code02.sh.e24507647
code02.sh.pe24507647
code02.sh.po24507647
code02.sh.o24507647 
list.txt
trimmed
adapters.fa
db                 
log.txt
fastqc                
salmon_result
```

### パート１終了。　次はRでの解析。
各ソフトのオプションについては以下を参照
<a href="https://qiita.com/gibukod/private/d9f3e656df3fdf9bde59">RNA-seq講習会＠水域保全学研究室　実践編補足資料</a>

**<a href="https://qiita.com/gibukod/private/dadb8615d34c76b96dbd">パート２へ移動</a>**

