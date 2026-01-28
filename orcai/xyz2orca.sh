#!/bin/bash

TEMPLATE_FILE="orca_temp"
orca='/opt/orca-6.1.1/orca'

filelist=""  

files=${filelist:-*.xyz}
for file in $files; do
    [[ -f $file ]] || continue
    # PYTHONWARNINGS="ignore" ase convert $file ${file%.*}.xyz
    sed "/\* XYZFILE/s|$| ${file%.*}.xyz|" "$TEMPLATE_FILE" > "${file%%.*}.freq.inp"
    # $orca ${file%.*}.freq.inp |tee ${file%.*}.freq.out
done

