#!/bin/bash
#BSUB -q mpi -n 24 -o %J.log -e %J.err -J job

# version 251025 by lgao1024@outlook.com
########################################## setup environments and functionals #########################################
module purge 
module load gcc/gcc-14.10 openmpi/4.1.4_gcc14
#export OMPI_MCA_osc=^ucx
ORCA_PATH=/share/home/syuan/lgao/apps/orca-6.1.0/
orca='/share/home/syuan/lgao/apps/orca-6.1.0/orca'
export PATH=$ORCA_PATH:$PATH
export LD_LIBRARY_PATH=$ORCA_PATH:$LD_LIBRARY_PATH

function mpi_orca() {
  local file=$1
  mkdir -p "${file%.*}" && mv -f "${file%%.*}".*.* "${file%%.*}"_*.* "${file%.*}/" 2>/dev/null
  ( cd "${file%.*}" && "$orca" "$file" > "${file%.*}.out" )
}

function mpi_orca_here() {
  local file=$1
  "$orca" "$file" > "${file%.*}.out"
}
################################################## start calculations  ################################################

filelist=""  

files=${filelist:-*.inp}
for file in $files; do
    [[ -f $file ]] || continue
    mpi_orca "${file}"
done
