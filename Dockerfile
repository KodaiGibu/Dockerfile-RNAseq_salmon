FROM ubuntu:20.04

RUN ln -sf /usr/share/zoneinfo/Asia/Tokyo /etc/localtime
RUN apt update -y && apt install -y fastqc python3 python-is-python3 python3-pip wget salmon trimmomatic vim
WORKDIR /bin
RUN wget https://github.com/shenwei356/seqkit/releases/download/v2.5.0/seqkit_linux_amd64.tar.gz && tar xvf seqkit_linux_amd64.tar.gz && rm seqkit_linux_amd64.tar.gz
RUN mkdir -p /home/app
WORKDIR /home/app
COPY ./code01.sh ./
