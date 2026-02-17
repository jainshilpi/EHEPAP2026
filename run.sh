#!/bin/bash

# Loop with different integer and float values
#for Nisgma in {1..5}; do

inputFilesList=inputfile
for Nsigma in $(seq 1 0.5 5); do
#for Nsigma in $(seq 5 1 10); do
    #    for pedwidth in $(seq 0.1 0.1 1.0); do ### seq startingVal step finalVal
    #        for pedwidth in $(seq 0.04 0 0.04); do
    pedwidth=0.04

    # Call your executable script with the current integer and float values
    echo "Value of Nsigma is "$Nsigma
    ./serc19_ecal_clustering "$Nsigma" "$pedwidth" "$inputFilesList"
    #done
done

