#!/bin/sh

debug_flag=0
clean_build=0

# check input arguments
for arg in "$@"
do
    if [ "$arg" == "-g" ] || [ "$arg" == "-debug" ]
    then
        debug_flag=1
    elif [ "$arg" == "-clean" ]
    then
        clean_build=1
    fi
done

# remove old files for clean build
if [ $clean_build == 1 ]
then
    rm -r build/*
fi

#build and install
cd build
if [ $debug_flag == 1 ]
then
    cmake -DCMAKE_BUILD_TYPE="Debug" -DCMAKE_INSTALL_PREFIX=$PWD/platform ..
else
    cmake -DCMAKE_BUILD_TYPE="Release" -DCMAKE_INSTALL_PREFIX=$PWD/platform ..
fi

make all
make install
