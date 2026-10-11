#!/bin/bash
mkdir -p results workdirs

for json in transformations/*.json; do
    name=$(basename "$json" .json)
    out="results/${name}_results.json"
    wdir="workdirs/$name"
    
    if [ -f "$out" ]; then
        echo "Skipping $name (already done)"
        continue
    fi
    
    echo "Running $name..."
    mkdir -p "$wdir"
    rm -rf "$wdir"/scratch_*
    
    taskset -c 0-3 openfe quickrun "$json" -o "$out" -d "$wdir"
done
