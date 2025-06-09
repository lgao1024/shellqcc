#!/bin/bash
#BSUB -q fat_384
#BSUB -o %J.log -e %J.err
#BSUB -n 80
#BSUB -J qec_task

# module purge 
# rm -f *.log *.err
# . /share/home/syuan/lgao/apps/intel2019u5.sh
# export PATH=$PATH:/share/home/syuan/lgao/apps/qe-7.3.1_intel2019u5/bin 
# export ESPRESSO_PSEUDO=/share/home/syuan/lgao/apps/PSEUDO/PBEsol
# export ESPRESSO_TMPDIR=/scratch/syuan/qe

# version 0.1-250609
READ_FILE="$1"
IN_FILE=$READ_FILE

run_qe () {

    # check_nbnd $IN_FILE
    calc_pw $IN_FILE
    # calc_pdos $IN_FILE
    # calc_co $IN_FILE
    # calc_esp $IN_FILE
    # calc_hubbard $IN_FILE
    # calc_bands $IN_FILE
    # calc_w90 $IN_FILE
    # calc_phonon $IN_FILE

}

check_nbnd () {
    
    IN_FILE="$1"
    DRY_TEMP=${IN_FILE%.*}.dry.tmp
    DRY_OUT=${IN_FILE%.*}.dry.out

    cp $IN_FILE $DRY_TEMP

    if grep -q "^ *electron_maxstep *= *" "$DRY_TEMP"; then
        sed -i "s/^ *electron_maxstep *= .*/  electron_maxstep = 0/" "$DRY_TEMP"
    else
        sed -i "/&ELECTRONS/a\  electron_maxstep = 0" "$DRY_TEMP"
    fi 

    mpirun pw.x -in $DRY_TEMP |tee $DRY_OUT


    suggested_band=$(awk -v ne="$nelec" 'BEGIN { printf "%d", int(ne/2*1.25) + 4 }')
    nks_num=$(grep -m1 -e 'number of Kohn-Sham states' "$DRY_OUT" | awk -F'=' '{print int($NF)}')
    band_num=$(( nks_num >= suggested_band ? nks_num : suggested_band ))

    echo "Number of electrons:      $nelec"
    echo "Suggested NBANDS:         $suggested_band"
    echo "Kohn-Sham states (nks):   $nks_num"
    echo "Final NBANDS (band_num):  $band_num"

    cp $IN_FILE $IN_FILE.bak && rm $DRY_TEMP

    if grep -q "^ *nbnd *= *" "$IN_FILE"; then
        sed -i "s/^ *nbnd *= .*/  nbnd = $band_num/" "$IN_FILE"
    else
        sed -i "/ntyp/a\  nbnd = $band_num" "$IN_FILE"
    fi 
}

calc_pw () {
    IN_FILE="$1"
    mpirun pw.x -in $IN_FILE |tee ${IN_FILE%.*}.out
}

calc_pdos () {

    IN_FILE="$1"

    if [[ ! -f ${IN_FILE%%.*}.dos.in ]]; then
        {
            echo "&DOS"
            grep -E "^\s*(prefix|outdir)\s*=" $IN_FILE
            echo "  degauss = 0.005"
            echo "/"
            echo "&PROJWFC"
            grep -E "^\s*(prefix|outdir)\s*=" $IN_FILE
            echo "  degauss = 0.005"
            echo "/"
        } > ${IN_FILE%%.*}.dos.in
    fi

    mpirun dos.x -in ${IN_FILE%%.*}.dos.in |tee ${IN_FILE%%.*}.dos.out
    
    mpirun projwfc.x -in ${IN_FILE%%.*}.dos.in |tee ${IN_FILE%%.*}.pdos.out

    for elem in $(grep "), wfc" ${IN_FILE%%.*}.pdos.out | awk -F'[()]' '{print $2}' | awk '{print $1}' | sort -u); do
        sumpdos.x *\($elem\)_* > ${IN_FILE%%.*}.pdos_$elem
    done

    mkdir -p ${IN_FILE%%.*}.pdos && mv ${IN_FILE%%.*}.pdos_* ${IN_FILE%%.*}.pdos/
}

calc_co () {

    IN_FILE="$1"
    hoco_num=$(grep -m1 -e 'number of electrons' "${IN_FILE%.*}.out" | awk '{print int($NF/2)}')
    co_min=$((hoco_num - 2))
    co_max=$((hoco_num + 3))

    echo "HOCO number:    $hoco_num"
    if [[ ! -f ${IN_FILE%%.*}.pp_co.in ]]; then
        {
            echo "&INPUTPP"
            grep -E "^\s*(prefix|outdir)\s*=" $IN_FILE
            echo "  filplot='${IN_FILE%%.*}.pp_co'"
            echo "  plot_num = 7"
            echo "  kpoint = 1"
            echo "  kband(1)= $co_min"
            echo "  kband(2)= $co_max"
            echo "  lsign = .true."
            echo "/"
            echo "&PLOT"
            echo "  iflag = 3"
            echo "  output_format = 6 "
            echo "  fileout = '.cube'"
            echo "/"
        } > ${IN_FILE%%.*}.pp_co.in
    fi

    mpirun pp.x -in ${IN_FILE%%.*}.pp_co.in |tee ${IN_FILE%%.*}.pp_co.out

    mkdir -p ${IN_FILE%%.*}.pp_co && mv ${IN_FILE%%.*}.pp_co_* ${IN_FILE%%.*}.pp_co/
}

calc_esp () {

    IN_FILE="$1"

    if [[ ! -f ${IN_FILE%%.*}.pp_esp.in ]]; then
        {
            echo "&INPUTPP"
            grep -E "^\s*(prefix|outdir)\s*=" $IN_FILE
            echo "  filplot='${IN_FILE%%.*}.pp_esp'"
            echo "  plot_num = 11"
            echo "/"
            echo "&PLOT"
            echo "  iflag = 3"
            echo "  output_format = 6"
            echo "  fileout = '${IN_FILE%%.*}.pp_esp.cube'"
            echo "  ! iflag = 1"
            echo "  ! output_format = 0"
            echo "  ! filplot(1)='${IN_FILE%%.*}.pp_esp'"
            echo "  ! fileout = '${IN_FILE%%.*}.pp_esp.gp'"
            echo "  ! e1(1) = 1.0, e1(2) = 0.0, e1(3) = 0"
            echo "  ! e2(1) = 0.0, e2(2) = 1.0, e2(3) = 0"
            echo "  ! e3(1) = 0.0, e3(2) = 0.0, e3(3) = 1.0"
            echo "  ! x0(1) = 0.5, x0(2) = 0.5, x0(3) = 0.5"
            echo "  nx = 100, ny = 100"
            echo "/"
        } > ${IN_FILE%%.*}.pp_esp.in
    fi

    mpirun pp.x -in ${IN_FILE%%.*}.pp_esp.in |tee ${IN_FILE%%.*}.pp_esp.out

}

calc_hubbard () {

    IN_FILE="$1"

    if ! grep -qi "HUBBARD" "$IN_FILE"; then
        echo "HUBBARD {ortho-atomic}" >> "$IN_FILE"
        echo "You need to add the Hubbard card like the following and re-run SCF:"
        echo "HUBBARD {ortho-atomic}"
        echo "U     Ti-3d  4"
        exit 1
    fi

    if [[ ! -f ${IN_FILE%%.*}.hp.in ]]; then
        {
            echo "&INPUTHP"
            grep -E "^\s*(prefix|outdir)\s*=" $IN_FILE
            echo "  nq1 = 1, nq2 = 1, nq3 = 1"
            echo "  conv_thr_chi = 1.0d-5"
            echo "  iverbosity = 2"
            echo "/"
        } > ${IN_FILE%%.*}.hp.in 
    fi

    mpirun pw.x -in $IN_FILE |tee ${IN_FILE%.*}.out
    mpirun hp.x -in ${IN_FILE%%.*}.hp.in |tee ${IN_FILE%%.*}.hp.out
}

calc_bands () {
    
    IN_FILE="$1"
    BANDS_IN=${IN_FILE%%.*}.bands.in

    if [[ ! -f ${IN_FILE%%.*}.band.in ]]; then
        {
            echo "&BANDS"
            grep -E "^\s*(prefix|outdir)\s*=" $IN_FILE
            echo "  filband ='${IN_FILE%%.*}.band'"
            echo "  lsym = .false."
            echo "/"
        } > ${IN_FILE%%.*}.band.in
    fi

    if [[ ! -f ${IN_FILE%%.*}.band_plot.in ]]; then
        Efermi=$(grep "highest occupied" ${IN_FILE%.*}.out | awk -F'[:]' '{print $2}' | awk '{print $1}')
        Eband_min=$(echo "$Efermi - 10" | bc)
        Eband_max=$(echo "$Efermi + 10" | bc)
        {
            echo "${IN_FILE%%.*}.band"
            echo "$Eband_min, $Eband_max"
            echo "${IN_FILE%%.*}.band_plot.xmgr"
            echo "${IN_FILE%%.*}.band_plot.ps"
            echo "$Efermi eV"
            echo "2 $Efermi eV"
        } > ${IN_FILE%%.*}.band_plot.in
    fi

    if [[ ! -f $BANDS_IN ]]; then

        awk '
            /Real form of k-point coordinates/ {flag=1; next}
            /^#/ || NF==0 {next}
            flag {
                if ($1 ~ /^#/) next
                count++
                x[count] = $1
                y[count] = $2
                z[count] = $3
                label[count] = $4
            }
            /END of FILE/ {flag=0}
            END {
                print count
                interval = int(100 / count)
                for (i = 1; i <= count; i++) {
                    printf "   %.10f     %.10f     %.10f     %3d   ! %s\n", x[i], y[i], z[i], interval, label[i]
                }
            }
        ' "${IN_FILE%%.*}.kpf" > band_kps.tmp

        awk -v kpath_file="band_kps.tmp" '
            BEGIN {
                inserting = 0
                skip = 0
            }
            /^\s*calculation\s*=/ {
                print "  calculation = '\''bands'\''"
                next
            }
            /^K_POINTS/ {
                inserting = 1
                if ($2 == "automatic") skip = 1
                print "K_POINTS {crystal_b}"
                # Read and print the k-point path from the external file
                while ((getline line < kpath_file) > 0) print line
                close(kpath_file)
                next
            }
            skip {
                skip = 0
                next
            }
            {
                print
            }
        ' "$IN_FILE" > "$BANDS_IN" && rm band_kps.tmp
    fi

    mpirun pw.x -in $BANDS_IN |tee ${BANDS_IN%.*}.out
    mpirun bands.x -in ${IN_FILE%%.*}.band.in |tee ${IN_FILE%%.*}.band.out
    plotband.x < ${IN_FILE%%.*}.band_plot.in |tee ${IN_FILE%%.*}.band_plot.out
    mkdir -p ${IN_FILE%%.*}.band_plot && mv ${IN_FILE%%.*}.band_plot.* ${IN_FILE%%.*}.band_plot/

}

calc_w90 () {

    IN_FILE=$1
    W90IN=${IN_FILE%%.*}.win

    if [[ ! -f ${IN_FILE%%.*}.open_grid.out ]]; then
        {
            echo "&INPUTPP"
            grep -E "^\s*(prefix|outdir)\s*=" $IN_FILE
            echo "/"
        } > ${IN_FILE%%.*}.open_grid.in

         mpirun open_grid.x -in ${IN_FILE%%.*}.open_grid.in |tee ${IN_FILE%%.*}.open_grid.out
    fi

    
    if [[ ! -f $W90IN ]]; then
        awk '
            /wannier90/ {flag=1; next}
            /OPEN_GRID    :/ {flag=0; next}
            flag && NF
        ' ${IN_FILE%%.*}.open_grid.out > w90_kpoints.tmp
        
        awk '
            /ATOMIC_POSITIONS/ {flag=1; next}
            /K_POINTS/ {flag=0; next}
            flag && NF
        ' $IN_FILE > w90_pos.tmp

        awk '
            $1 == "A" {
                a = $3
            }
            $1 == "CELL_PARAMETERS" {
                flag = 1
                next
            }
            $1 == "ATOMIC_SPECIES" {
                flag = 0
            }
            flag && NF {
                for (i = 1; i <= NF; i++) {
                    v[i] = $i * a
                }
                printf "   %.7f    %.7f    %.7f\n", v[1], v[2], v[3]
            }
        ' "$IN_FILE" > w90_cell.tmp

        awk '
            /Real form of k-point coordinates/ { flag = 1; next }
            /^#/ || NF == 0 { next }
            flag && $1 !~ /^#/ {
                labels[++count] = $NF
                coords[count] = sprintf("%.10f %.10f %.10f", $1, $2, $3)
            }
            /END of FILE/ { flag = 0 }
            END {
                for (i = 1; i < count; i++) {
                    printf "%-6s  %s    %-6s  %s\n", labels[i], coords[i], labels[i+1], coords[i+1]
                }
            }
        ' ${W90IN%.*}.kpf > w90_kpath.tmp

        mp_grid=$(awk '
                     tolower($0) ~ /^k_points/ && tolower($2) == "automatic" { getline; print $1, $2, $3 }
                ' "$IN_FILE")

        nbnd=$(grep -E "^\s*nbnd\s*=" "$IN_FILE" | awk -F'[=]' '{print int($NF)}')
        Efermi=$(grep "highest occupied" ${IN_FILE%.*}.out | awk -F'[:]' '{print $2}' | awk '{print $1}')
        Eband_min=$(echo "$Efermi - 10" | bc)
        Eband_max=$(echo "$Efermi + 10" | bc)
        Eband_froz_max=$(echo "$Efermi - 6" | bc)


        {
            echo "num_bands = $nbnd"
            echo "num_wann  = 20"
            # echo "exclude_bands = 1-12"
            echo "dis_win_min = $Eband_min"
            echo "dis_win_max = $Eband_max"
            echo "dis_froz_max = $Eband_froz_max"
            # echo "conv_window = 3"
            # echo "conv_tol = 1.0d-10"
            # echo "dis_num_iter = 200"
            # echo "num_iter = 100"
            echo "bands_plot = .true."
            echo "bands_num_points = 100"
            # echo "guiding_centres = .true."
            echo "auto_projections = .true."
            # echo ""
            # echo "begin projections"
            # echo "!random"
            # echo "O:px;py;pz"
            # echo "Ti:dz2;dx2-y2;dxz;dyz;dxy"
            # echo "end projections"
            echo ""
            echo "Begin Kpoint_Path"
            cat w90_kpath.tmp
            echo "End Kpoint_Path"
            echo ""
            echo "begin unit_cell_cart"
            echo "ang"
            cat w90_cell.tmp
            echo "end unit_cell_cart"
            echo ""
            echo "begin atoms_frac"
            cat w90_pos.tmp
            echo "end atoms_frac"
            echo ""
            echo "mp_grid: $mp_grid"
            echo "begin kpoints"
            cat w90_kpoints.tmp
            echo "end kpoints"
        } > $W90IN && rm w90_*.tmp
    fi

    if [[ ! -f ${IN_FILE%%.*}.pw2wan.in ]]; then
        {
            echo "&INPUTPP"
            grep -E "^\s*prefix\s*=" "$IN_FILE" | sed -E "s/(prefix\s*=\s*')([^']+)(')/\1\2_open\3/"
            grep -E "^\s*outdir\s*=" "$IN_FILE"
            echo "  seedname = '${W90IN%.*}'"
            # echo "  atom_proj = .true."
            echo "  scdm_proj = .true."
            echo "  scdm_entanglement = 'erfc'" 
            echo "  scdm_mu = 10"
            echo "  scdm_sigma = 4"
            echo "  write_unk = .true."
            echo "/"
        } > "${IN_FILE%%.*}.pw2wan.in"
    fi

    wannier90.x -pp $W90IN
    mpirun pw2wannier90.x -in ${IN_FILE%%.*}.pw2wan.in |tee ${IN_FILE%%.*}.pw2wan.out
    wannier90.x $W90IN
    mkdir -p ${IN_FILE%%.*}.w90 && mv UNK*.1 ${IN_FILE%%.*}.w90/
    gnuplot -e "set terminal pdfcairo; set output '${W90IN%.*}_band.pdf'" ${W90IN%.*}_band.gnu 
}

calc_phonon () {

    IN_FILE="$1"
    PH_BASENAME=${IN_FILE%%.*}
    mp_grid=$(awk '
                tolower($0) ~ /^k_points/ && tolower($2) == "automatic" { getline; print $1, $2, $3 }
            ' "$IN_FILE")
    read nk1 nk2 nk3 <<< "$mp_grid"

    if [[ ! -f $PH_BASENAME.ph.in ]]; then
        {
            echo "$PH_BASENAME"
            echo "&INPUTPH"
            grep -E "^\s*(prefix|outdir)\s*=" $IN_FILE
            echo "  dftd3_hess='automatic.hess'"
            echo "  fildyn='$PH_BASENAME.dyn'"
            echo "  tr2_ph = 1d-14"
            echo "  ! recover = .true."
            echo "  ldisp=.true."
            echo "  nq1=$nk1, nq2=$nk2, nq3=$nk3"
            echo "/"
        } > $PH_BASENAME.ph.in
    fi

    if [[ ! -f $PH_BASENAME.q2r.in ]]; then
        {
            echo "&INPUT"
            echo "  zasr='crystal'"
            echo "  fildyn='$PH_BASENAME.dyn'"
            echo "  flfrc='$PH_BASENAME.fc'"
            echo "/"
        } > $PH_BASENAME.q2r.in
    fi

    if [[ ! -f $PH_BASENAME.matdyn.in ]]; then
        {
            echo "&INPUT"
            echo "  asr='crystal'"
            echo "  flfrc='$PH_BASENAME.fc'"
            echo "  flfrq='$PH_BASENAME.freq'"
            echo "  flvec='$PH_BASENAME.modes'"
            echo "  dos   = .true."
            echo "  fldos = '$PH_BASENAME.dos'"
            echo "  nk1 = $nk1, nk2 = $nk2, nk3 = $nk3"
            echo "  ! q_in_band_form=.true."
            echo "/"
            echo "! 1"
            echo "! 0 0 0 1"
            echo "! 0 0 0 0"
        } > $PH_BASENAME.matdyn.in
    fi

    mpirun ph.x -in $PH_BASENAME.ph.in |tee $PH_BASENAME.ph.out
    mpirun q2r.x -in $PH_BASENAME.q2r.in |tee $PH_BASENAME.q2r.out
    mpirun matdyn.x -in $PH_BASENAME.matdyn.in |tee $PH_BASENAME.matdyn.out

}

run_qe
