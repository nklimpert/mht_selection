#!/bin/bash

dirs=$1

while read line; do
    mkdir $line
done < $dirs