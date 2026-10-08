#!/bin/bash

# 1. Definir carpetas base (sin '/' al inicio)
TRANS_DIR="transformations"
results="results"
workdirs="workdirs"

# 2. Crear carpetas principales
mkdir -p "$results" "$workdirs"

# 3. Iterar sobre las transformaciones
for json_file in "$TRANS_DIR"/*.json; do
    transname=$(basename "$json_file" .json)
    
    # Rutas relativas explícitas con ./
    out_json="./${results}/${transname}_results.json"
    target_workdir="./${workdirs}/${transname}"
    
    echo "==============================================================================="
    echo "               Running system: $transname"
    echo "==============================================================================="
    
    # Crear la subcarpeta de trabajo específica antes de ejecutar
    mkdir -p "$target_workdir"
    
    # Ejecutar OpenFE quickrun
    taskset -c 0-3 openfe quickrun "$json_file" \
        -o "$out_json" \
        -d "$target_workdir"
        
    echo "Results: $out_json"
    echo "Workdir: $target_workdir"
    echo "-----------------------------------------------------------------"
done
