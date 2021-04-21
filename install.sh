#!/bin/bash

set -e

#Install LRez
cd LRez
./install.sh

#Install LEVIATHAN
cd ../
make -j"$( nproc )" PREFIX="$( pwd )"
