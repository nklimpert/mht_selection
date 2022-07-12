#!/bin/bash

find paml/*/*/*.csv -exec cat {} \; >> combined_results.csv

python scripts/paml_stats.py combined_results.csv