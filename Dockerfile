############################################################
# Dockerfile to build cloneHD workflow container
# Based on Ubuntu
############################################################

FROM ubuntu

# File Author / Maintainer
MAINTAINER Ignacio Vazquez-Garcia <ivg@sanger.ac.uk>

# Install software
RUN apt-get update && apt-get install -y gfortran \
make gcc build-essential wget libgsl2 gsl-bin libgsl-dev libboost-all-dev \
libblas-dev liblapack-dev git perl python-pip gzip

WORKDIR /opt

# Install python modules
RUN pip install PyVCF

# Compile 
RUN make -f Makefile
# Install cloneHD
RUN git clone https://github.com/ivazquez/cloneHD.git && cd cloneHD && git checkout smchet
RUN cd cloneHD/src && mkdir ../build && make -f Makefile.farm

# Install smchet-challenge workflow
RUN git clone https://github.com/ivazquez/smchet-challenge.git && git checkout develop && cd smchet-challenge && make -f Makefile

# Copy scripts to `WORKDIR`
COPY smchet_workflow.sh *.pl *.py run_metrics ./