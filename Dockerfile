# Get the base Ubuntu image from Docker Hub
FROM ubuntu:latest

# Specify the working directory
WORKDIR /app

# Install dependencies on the base image
RUN apt-get -y update && apt-get install -y \
    g++ \
    libboost-all-dev \
    cmake

# Copy the current folder which contains C++ source code to the Docker image under /usr/src
COPY . /app

RUN rm -rf build

# Use Clang to compile the Test.cpp source file
RUN chmod +x build.sh && ./build.sh

# go back to app directory so that output files go there
RUN cd /app

# Run the output program from the previous step
ENTRYPOINT ["./build/pepsirf"]