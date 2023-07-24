# RNAseq_salmon
Dockerfile for RNA-seq with salmon

## Dokerfile:RNAseq_salmon
このDockerfileは、トランスクリプトームデータをリファレンスとしたRNA-seq解析用パイプラインに必要なソフトを格納したDockerfileである。
以下のソフトがインストールされている。
```sh:make image 
seqkit
Trimmomatic
fastqc
salmon
```
## imageの作成
リポジトリよりDockerfile、code01.shをダウンロードし、同じディレクトリに置く。


そして以下のコマンドを実行。
```sh:make image 
#docker imageの作成。rnaseq_salmonは自由表記。
docker image build -f Dockerfile -t rnaseq_salmon .

#imageが作成できたか確認
docker images
#問題なくイメージが作成されたら、以下のように表示される
REPOSITORY                  TAG             IMAGE ID       CREATED          SIZE
rnaseq_salmon               latest          66cb8ab17d4d   43 minutes ago   1GB
```

## 解析の準備

マウントするディレクトリに解析するfastqファイルおよびadapters.faを置く。

imageからdocker コンテナを立ち上げる。
以下ではsalmonという名前を付けたdocker コンテナを立ち上げている
以下では、home下にrnaseq_salmonというディレクトリを作成し、コンテナ内のmntとマウントしている。
```sh:dockerコンテナの立ち上げ
#コンテナの立ち上げ
docker run -itd -v /home/rnaseq_salmon/:/mnt --name salmon rnaseq_salmon /bin/bash

#コンテナに入る
docker exec -it salmon /bin/bash
```
以下のコマンドを実行してmntに移動し、fastqファイルおよびadapters.faがあるかを確認する。appにあるcode01.shをmntへコピーする。
```sh:dockerコンテナの立ち上げ
#mntへ移動
cd /mnt/

#mntの中身を確認
ls

#code01.shをコピー
cp /home/app/code01.sh ./
```

### リファレンスRNAデータのダウンロードおよびindexの作成
リファレンスとなるトランスクリプトームデータをダウンロードし、salmon用のindexを作成する。
以下は例としてオニヒトデのトランスクリプトームデータをダウンロードする。
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

### indexのパスcode01.shに書き込む
pwd等のコマンドでindexのディレクトリを確認し、code01.shの`db="/home/gibukodai/216/Ap/db/Apla"`のpathを書き換える。adapters.faを編集し、シーケンスに利用したアダプター配列の情報を記載する。
pwdでは、indexの置かれたディレクトリまでは出ないので注意。
```ディレクトリ内のファイルの確認
#pathの確認。
pwd

#code01.shを編集
cd /mnt/
vim code01.sh

#adapters.faを編集
vim adapters.fa
```
adapters.faはアダプター配列をfastaファイルとして保存たものである。
```adapters.fa
>Read_1
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
>Read_2
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
```
## 解析の実行
以下のようにスクリプトを実行し、解析が終わるまで待つ。

```sh:スクリプトの実行
#以下のコードによりバックグラウンドでジョブを実行
nohup bash code01.sh &
#実行を確認
ps -u



#スクリプトに許可がない場合は実行できないので、以下のようにして許可を与える(permission deniedと表示される)
chmod +x code01.sh

```

<br>

## 出力ファイル
実行が終わると以下のようにファイルが出力されている。
Sum1.txtにrawdataのリード数、trimmedにトリミング済みfastq、fastqcにクオリティの記されたhtmlファイル、salmon_resultにリードカウントの結果が格納されている。
次の解析にはsalmon_resultを用いる。


