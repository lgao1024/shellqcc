#!/bin/bash
#BSUB -q serial -n 1 -o %J.log -e %J.err -J zeopp

# version 251025 by lgao1024@outlook.com
module purge
module load gcc/gcc-14.10 openmpi/4.1.4_gcc14
ZEO_BIN='/share/home/syuan/lgao/apps/zeo++-0.3/network'

function run_zeopp() {
    local file=$1
    # $ZEO_BIN -r rad.rad -mass mass.mass -psd 1.86 1.86 3000 ${file}
    # $ZEO_BIN -ha -psd 1.86 1.86 3000 ${file}
    # $ZEO_BIN -gridG ${file} > grid.out

    $ZEO_BIN -ha -resex -volpo 1.86 1.86 30000 -sa 1.86 1.86 6000 -psd 1.86 1.86 3000 ${file}
    echo "" >> summary.out
    echo "#" $file >> summary.out
    grep -e 'PROBE_OCCUPIABLE___RESULT' ${file%.*}.volpo | awk '{print "Ads(cm^-3/g STP): " $(NF-3)*0.8064/$3/28.1*22.4*1000}' >> summary.out
    grep @ ${file%.*}.sa | awk '{print "SA(cm^-2/g): " $(NF-6)}' >> summary.out
}



filelist=""  

files=${filelist:-*.cif}
for file in $files; do
    [[ -f $file ]] || continue
    run_zeopp $file
done
