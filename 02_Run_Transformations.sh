#!/bin/bash

mkdir -p results workdirs

for json in transformations/*.json; do
    name=$(basename "$json" .json)
    out="results/${name}_results.json"
    wdir="workdirs/$name"
    
    # Check if already completed (Resume logic)
    if [ -f "$out" ]; then
        echo "Skipping: $name (already completed)"
        continue
    fi
    
    echo "========================================================="
    echo "   Running system: $name"
    echo "========================================================="
    
    # Delete the entire folder to remove ALL clutter from failed attempts, then recreate
    rm -rf "$wdir"
    mkdir -p "$wdir"
    
    # Run OpenFE
    taskset -c 0-3 openfe quickrun "$json" -o "$out" -d "$wdir"
        
    echo "Results: $out"
    echo "Workdir: $wdir"
    echo "---------------------------------------------------------"
done
