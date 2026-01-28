#!/bin/bash

TEMPLATE_FILE="orca_temp_xtb"
orca='/opt/orca-6.1.1/orca'

filelist=""  

files=${filelist:-*.mol}
for file in $files; do
    [[ -f $file ]] || continue
    PYTHONWARNINGS="ignore" ase convert $file ${file%.*}.xyz
    sed "/\* XYZFILE/s|$| ${file%.*}.xyz|" "$TEMPLATE_FILE" > "${file%.*}.xtb.inp"
    $orca ${file%.*}.xtb.inp |tee ${file%.*}.xtb.out
done

