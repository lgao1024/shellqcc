#!/bin/bash
#BSUB -q mpi -n 24 -o %J.log -e %J.err -J job

. /share/home/syuan/lgao/apps/intel2019u5.sh
export PATH=$PATH:/share/home/syuan/lgao/apps/qe-7.5/bin 
export ESPRESSO_PSEUDO=/share/home/syuan/lgao/apps/PSEUDO/upf_PBEsol
export ESPRESSO_TMPDIR=/scratch/syuan/qe


qec.sh TiErBDCH.scf.in

