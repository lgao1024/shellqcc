#!/bin/bash
#BSUB -q mpi -n 24 -o %J.log -e %J.err -J job

# version 251025 by lgao1024@outlook.com
########################################## setup environments and functionals #########################################
module purge 
module load gcc/gcc-14.10 openmpi/4.1.4_gcc14
#export OMPI_MCA_osc=^ucx
source /share/home/syuan/lgao/apps/cp2k-2025.2/tools/toolchain/install/setup
CP2K_PATH=/share/home/syuan/lgao/apps/cp2k-2025.2/exe/local
export PATH=$CP2K_PATH:$PATH

function mpi_cp2k_here() {
    local file=$1
    mpirun cp2k.psmp -in "$file" > "${file%.*}.out"
}
################################################## start calculations  ################################################



filelist=""  

files=${filelist:-*.inp}
for file in $files; do
    [[ -f $file ]] || continue
    mpi_cp2k_here "${file}"
done