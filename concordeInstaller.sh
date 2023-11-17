#! /bin/bash

sudo ls     #useful to ask the sudo password immediatly rather than at the end risking to losing the data
wget https://www.math.uwaterloo.ca/tsp/concorde/downloads/codes/src/co031219.tgz
gunzip co031219.tgz
tar xvf co031219.tar
rm ./co031219.tar
mkdir concordeBuild && cd concordeBuild
CFLAGS="-g -O3 -march=native -mtune=native" ../concorde/configure
make
cd ..
sudo mv ./concordeBuild /opt/concorde
sudo rm -r ./concorde
echo ""
echo "ALL DONE"

