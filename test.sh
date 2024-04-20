#!/usr/bin/env bash

PROC=4
FILENAME="$1"
STEPS="$2"

if [ "$#" -ne 2 ]; then
    echo "Illegal number of parameters"
    exit 1
fi

mpic++ --prefix /usr/local/share/OpenMPI -o life life.cpp
mpirun --prefix /usr/local/share/OpenMPI --use-hwthread-cpus -np "$PROC" life "$FILENAME" "$STEPS"

