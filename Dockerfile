FROM ubuntu:latest

LABEL maintainer="Sami Teeny <samiteeny@gmail.com>"

RUN apt-get update && apt-get install -y python python3

WORKDIR /pathogenicity
ADD ./pathogenicity .
WORKDIR /minpath
ADD ./minpath .
WORKDIR /usr/local/bin
ADD ./diamond .
ADD ./humiid .

ENTRYPOINT ["humiid"]
CMD ["help"]
