
FROM ubuntu:18.04

LABEL org.opencontainers.image.authors="Camilla FERRARI"

ADD requirements.txt requirements.txt

# Installing python 3.6 and pip3
RUN apt-get update
RUN apt-get install -y python3.6 python3-pip 
RUN pip3 install --upgrade pip

# Installing dependencies
RUN pip3 install -r requirements.txt
