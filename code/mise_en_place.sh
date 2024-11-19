#!/bin/sh

# mise_en_place.sh
# Benjamin D. Peterson

# This is the script that pulls down the scripts
# from Github that are needed to run this workflow.

cd ~/BLiMMP/code
git clone https://github.com/petersonben50/HomeBio.git


# Get updates
cd ~/BLiMMP/code/HomeBio
#git fetch
#git diff ...origin
git pull
rm -rf $CONDA_PREFIX/lib/python3.9/site-packages/homebio_functions
cp -r homebio_functions $CONDA_PREFIX/lib/python3.9/site-packages

