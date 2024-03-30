FROM ubuntu:22.04

# set the working directory in the container
WORKDIR /code

# install app dependencies
RUN apt-get update && apt-get install -y python3 python3-pip
RUN apt install python-is-python3
RUN apt-get install gcc
RUN apt-get install -y openmpi-bin libopenmpi-dev 

# copy the dependencies file to the working directory
COPY requirements.txt .

# install dependencies
RUN pip install -r requirements.txt

# copy the content of the local src directory to the working directory
RUN mkdir dsmpy
COPY dsmpy/ dsmpy/

# build fortran lib
COPY build.py .
RUN pwd
RUN ls
RUN python3 build.py

# copy the content of the local src directory to the working directory
COPY tests/ .
