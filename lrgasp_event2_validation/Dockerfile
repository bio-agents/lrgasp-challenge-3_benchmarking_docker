# Use an official Python runtime as a parent image
FROM python:3.7-slim

# Set the working directory to /app
WORKDIR /app

# Install any needed packages specified in requirements.txt
RUN apt-get update  && apt-get -y install procps && rm -rf /var/lib/apt/lists/*
RUN /usr/local/bin/python -m pip install --upgrade pip
RUN pip install jsonschema
RUN pip install numpy
RUN pip install fasta-reader
RUN pip install gtfparse
RUN pip install validators

# Copy the current directory contents into the container at /app
COPY . /app


