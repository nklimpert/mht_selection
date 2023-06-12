#!/bin/bash

for dir in relax/*/; do
    cd $dir;
    for subdir in */; do
        cd $subdir
        hyphymp relax_run.bf
        cd ..
    done
    cd ../../
done