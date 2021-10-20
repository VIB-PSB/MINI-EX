
FROM python:3.6

LABEL org.opencontainers.image.authors="Camilla FERRARI"

ADD requirements.txt requirements.txt

RUN pip install -r ./requirements.txt
