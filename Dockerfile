FROM ubuntu

RUN apt-get update && apt-get install -y gfortran build-essential \
make gcc build-essential wget libgsl0ldbl gsl-bin libgsl0-dev git \
libblas-dev liblapack-dev python-pip

WORKDIR /opt

RUN git clone -b smchet https://github.com/ivazquez/cloneHD.git
RUN cd cloneHD/src && mkdir ../build && make -f Makefile.farm

RUN git clone -b smchet https://github.com/ivazquez/cloneHD-tools.git
RUN cd cloneHD-tools && python setup.py install
