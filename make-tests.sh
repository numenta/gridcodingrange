#!/bin/bash

set -e

cd "$(dirname "$0")"

debug=false

outbin="run-tests"

cmd="g++ -o $outbin ./src/test/main.cpp ./src/test/GridUniquenessTest.cpp ./src/GridUniqueness.cpp ./src/external/gtest/src/gtest-all.cc -I./src -I./src/external -I./src/external/gtest -lpthread"

if [ "$debug" = true ] ; then
    cmd="$cmd -g -D NTA_ASSERTIONS_ON"
else
    cmd="$cmd -O3"
fi

eval $cmd

echo "To run unit tests, execute: ./$outbin"
