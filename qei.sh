#!/bin/bash

# version 0.1-260109-6
ASE_PATH="ase"
CIF2CELL_PATH="cif2cell"
CIF2CELL_IN="cif2cell.in"
SCF_K_RESOLUTION="0.35"
PP_FILE_NAME="-pbesol-std-oncv.UPF"
DEGAUSS="0.001"
PSEUDO_DIR="./pseudo"
OUT_DIR="./tmp"
PSEUDO_FUNC="pbesol"
PSEUDO_PATH=/opt/_pseudo/PBEsol
export VASP_PAWLIB=/opt/_pseudo/paw64_PBE/
export VASP_PP_PRIORITY="_d,_pv,_sv,_h,_s"


cif2qe() {

    CIF_FILE="$1"
    CALC_TYPE="$2"

    $CIF2CELL_PATH "$CIF_FILE" -p pwscf -o $CIF2CELL_IN \
        --setup-all --no-reduce --print-digits=10 \
        --k-resolution=$SCF_K_RESOLUTION \
        --pwscf-pseudostring=$PP_FILE_NAME 

    # &CONTROL block
    awk -v calc_type="$CALC_TYPE" \
        -v prefix="${CIF_FILE%.*}" \
        -v pseudo_dir="$PSEUDO_DIR" \
        -v out_dir="$OUT_DIR" \
        'BEGIN {
            print "&CONTROL"
            print "  calculation = '\''" calc_type "'\''"
            print "  prefix = '\''" prefix "'\''"
            print "  !verbosity = '\''high'\''"
            print "  !pseudo_dir = '\''" pseudo_dir "'\''"
            print "  outdir = '\''" out_dir "'\''"
            print "  !disk_io = '\''low'\''"
            print "  !restart_mode = '\''restart'\''"
            if (calc_type == "relax" || calc_type == "vc-relax") {
                print "  nstep = 200"
                print "  tprnfor = .true."
                print "  tstress = .true."
                print "  forc_conv_thr = 1.0d-3"
                print "  etot_conv_thr = 1.0d-4"
            }
            print "/"
        }
        NR > 9 { print }' "$CIF2CELL_IN" > temp && mv temp "$CIF2CELL_IN"

    # &SYSTEM block
    read -p "Set charge and magnetization [0 0]: " CHARGE MAG; CHARGE=${CHARGE:-0}; MAG=${MAG:-0}
    read -p "Assigin the DFT functional: " DFTINPUT VDWCORR; DFTINPUT=${DFTINPUT:-n}; VDWCORR=${VDWCORR:-n}
    read -p "Set ecutwfc for the pseudo [40 Ry]: " ECUT; ECUT=${ECUT:-40}

    awk -v charge="$CHARGE" \
        -v magnetization="$MAG" \
        -v dftinput="$DFTINPUT" \
        -v vdwcorr="$VDWCORR" \
        -v ecut="$ECUT" \
        'BEGIN { in_system = 0; done = 0 }
        /^&SYSTEM/ { in_system = 1 }
        in_system && /^\/$/ && !done {
            print "  ecutwfc = " ecut
            print "  ecutrho = " (ecut * 4)
            if (charge != "0") print "  tot_charge = " charge
            if (magnetization != "0") {
                print "  nspin = 2"
                print "  tot_magnetization = " magnetization
                print "  occupations = '\''smearing'\''"
                print "  smearing = '\''mv'\''"
                print "  degauss = 0.001"
                print "  !lspinorb = .true."
                print "  !noncolin = .true."
            }

            if (dftinput != "n") {
                print "  input_dft = '\''" (dftinput == "hsesol" ? "sla+pw+hse+psc" : dftinput) "'\''"
            }
            if (vdwcorr == "d3" || vdwcorr == "d3bj") {
                print "  vdw_corr = '\''dft-d3'\''"
                print "  dftd3_version = " (vdwcorr == "d3" ? 3 : 4)
            }

            done = 1
        }
        { print }
        ' "$CIF2CELL_IN" > temp && mv temp "$CIF2CELL_IN"

    # &ELECTRONS, &IONS, and &CELL blocks
    awk -v calc_type="$CALC_TYPE" '
        /^&SYSTEM/ { in_system = 1 }
        /^\/$/ {
            print
            if (in_system) {
                print "&ELECTRONS"
                print "  mixing_beta = 0.7"
                print "  conv_thr = 1.0d-6"
                print "  electron_maxstep = 200"
                print "  !scf_must_converge = .false."
                print "  !diagonalization = '\''cg'\''"
                print "  !startingwfc = '\''random'\''"
                print "/"
                if (calc_type == "relax" || calc_type == "vc-relax") {
                    print "&IONS"
                    print "  ion_dynamics = '\''bfgs'\''"
                    print "/"
                }
                if (calc_type == "vc-relax") {
                    print "&CELL"
                    print "  cell_dynamics = '\''bfgs'\''"
                    print "  press_conv_thr = 0.5d0"
                    print "  cell_dofree = '\''all'\''"
                    print "/"
                }
                in_system = 0
                next
            }
            next
        }
        { print }
        ' "$CIF2CELL_IN" > temp && mv temp "$CIF2CELL_IN"

    # ATOMIC_SPECIES block
    read -p "Check target pseudopotential type [US/paw/nc]: " PSEUDO_TYPE; PSEUDO_TYPE=${PSEUDO_TYPE:-us}
    read -p "Copy pseudo files to ${PSEUDO_DIR} ? [y/N] " CP_PSEUDO; CP_PSEUDO=${CP_PSEUDO:-n}
    
    mapfile -t ELEM_LIST < <(grep -o '\S*\.UPF' "$CIF2CELL_IN" | awk -F'-' '{print $1}')
    for ELEM in "${ELEM_LIST[@]}"; do
        files=($(ls $PSEUDO_PATH/${ELEM}-*${PSEUDO_FUNC}-*${PSEUDO_TYPE}*.UPF 2>/dev/null))

        if [ ${#files[@]} -eq 0 ]; then
            echo "No pseudopotential files matching current criteria"
            exit 1
        elif [ ${#files[@]} -eq 1 ]; then
            PSEUDO_NEW=$(basename "${files[0]}")
        else
            echo "Setup pseudo file for $ELEM :"
            for i in "${!files[@]}"; do
                echo "   $i: $(basename "${files[$i]}")"
            done
            read -p ">[0] " PSEUDO_CHOICE
            PSEUDO_CHOICE=$([[ "$PSEUDO_CHOICE" =~ ^[0-9]+$ ]] && ((PSEUDO_CHOICE >= 0 && PSEUDO_CHOICE < ${#files[@]})) && echo "$PSEUDO_CHOICE" || echo 0)
            PSEUDO_NEW=$(basename "${files[$PSEUDO_CHOICE]}")
        fi

        sed -i "s/$(grep -o '\S*\.UPF' "$CIF2CELL_IN" | grep "${ELEM}-")/${PSEUDO_NEW}/g" "$CIF2CELL_IN"
        
        if [[ "$CP_PSEUDO" =~ ^[Yy]$ ]]; then
            mkdir -p "$PSEUDO_DIR"
            cp "$PSEUDO_PATH/$PSEUDO_NEW" "$PSEUDO_DIR"
        fi
    done

    # K_POINTS block
    awk '
        /^K_POINTS/ {
            getline nextline
            if (nextline ~ /^1[[:space:]]+1[[:space:]]+1[[:space:]]+0[[:space:]]+0[[:space:]]+0$/) {
                print "K_POINTS gamma"
                next
            } else {
                print $0
                print nextline
                next
            }
        }
        { print }
        ' "$CIF2CELL_IN" > temp && mv temp "$CIF2CELL_IN"
    
    mv $CIF2CELL_IN ${CIF_FILE%.*}.$CALC_TYPE.in
}  

cif2vasp () {
    
    CIF_FILE="$1"

    mkdir -p ${CIF_FILE%.*} && cd ${CIF_FILE%.*}
    $CIF2CELL_PATH "../$CIF_FILE" -p vasp --vasp-format=5 \
        --setup-all --no-reduce --print-digits=10 \
        --k-resolution=$SCF_K_RESOLUTION  
    NBANDS=$(grep NBANDS INCAR | awk '{print $NF}')
    
    awk -v nbands="$NBANDS" \
        'BEGIN {
            print "&CONTROL"
            print "ISTART = 1"
            print "ISPIN = 1"
            print "LREAL = Auto"
            print "# PREC = Accurate"
            print "# ADDGRID = .TRUE."
            print "# LH5 = .TRUE."
            print "# LWAVE = .FALSE."
            print "# LCHARG = .FALSE."
            print "# LVTOT  = .TRUE."
            print "NCORE = 4"
            print ""
            print "&ELECTRONS"
            print "ENCUT = 450"
            print "GGA = PS"
            print "IVDW = 12"
            print "NELMIN = 4"
            print "NELM = 200"
            print "EDIFF = 1E-05"
            print "# NBANDS = " nbands
            print "ISMEAR = 0"
            print "SIGMA = 0.05"
            print "NEDOS = 2001"
            print ""
            print "&IONS"
            print "NSW = 200"
            print "ISIF = 2                    (Stress/relaxation: 2-Ions, 3-Shape/Ions/V, 4-Shape/Ions)"
            print "IBRION = 1                  (Algorithm: 0-MD, 1-Quasi-New, 2-CG)"
            print "POTIM = 0.5"
            print "EDIFFG = -0.01"
        }' > INCAR && cd ..
}

cif2cp2k () {

    CIF_FILE="$1"
    cp2k_CALC_TYPE="GEO_OPT"

    $CIF2CELL_PATH "$CIF_FILE" -p cp2k -o $CIF2CELL_IN --no-reduce --print-digits=10

    # &GLOBAL and &FORCE_EVAL blocks
    awk -v proj="${CIF_FILE%.*}.$cp2k_CALC_TYPE" \
        'BEGIN {
            in_coord = 0
        }

        NR == 10 {
            print "&GLOBAL"
            print "  PROJECT " proj
            print "  PRINT_LEVEL LOW"
            print "  RUN_TYPE GEO_OPT"
            print "&END GLOBAL"
            print ""
            print "&FORCE_EVAL"
            print "  METHOD Quickstep"
            print "  &SUBSYS"
        }

        /&COORD/     { in_coord = 1 }
        /&END COORD/ { in_coord = 0 }

        NR > 9 {
            if (in_coord && NF == 4 && $1 ~ /^[A-Za-z]+$/) {
            elem = (length($1) == 1) ? " " $1 : $1
            printf("     %-2s %20.15f %20.15f %20.15f\n", elem, $2, $3, $4)
            if (!(elems_seen[$1]++)) elems[++n] = $1
            } else {
            print "    " $0
            }
        }

        END {
            print  ""
            for (i = 1; i <= n; i++) {

                printf "    &KIND %s\n", elems[i]
                printf "      ELEMENT %s\n", elems[i]
                print  "      BASIS_SET DZVP-MOLOPT-SR-GTH"
                print  "      POTENTIAL GTH-PBE"
                print  "    &END KIND"
            }
            print "  &END SUBSYS"
            print ""
            print "  &DFT"
            print "    BASIS_SET_FILE_NAME  BASIS_MOLOPT"
            print "    POTENTIAL_FILE_NAME  GTH_POTENTIALS"
            print "    CHARGE    0"
            print "    MULTIPLICITY    1"
            print "    &QS"
            print "      EPS_DEFAULT 1.0E-12"
            print "    &END QS"
            print "    &XC"
            print "      &XC_FUNCTIONAL"
            print "        &PBE"
            print "          PARAMETRIZATION PBESOL"
            print "        &END PBE"
            print "      &END XC_FUNCTIONAL"
            print "    &END XC"
            print "    &MGRID"
            print "      CUTOFF  400"
            print "      REL_CUTOFF  60"
            print "    &END MGRID"
            print "    &SCF"
            print "      MAX_SCF 200"
            print "      EPS_SCF 1.0E-06"
            print "#     SCF_GUESS RESTART"
            print "#     IGNORE_CONVERGENCE_FAILURE"
            print "      &PRINT"
            print "        &RESTART"
            print "          BACKUP_COPIES 0"
            print "        &END RESTART"
            print "      &END PRINT"
            print "      &DIAGONALIZATION"
            print "        ALGORITHM STANDARD"
            print "      &END DIAGONALIZATION"
            print "      &MIXING"
            print "        METHOD BROYDEN_MIXING"
            print "        ALPHA 0.4"
            print "        NBROYDEN 8"
            print "      &END MIXING"
            print "    &END SCF"
            print "  &END DFT"
            print "&END FORCE_EVAL"
            print ""
            print "&MOTION"
            print "  &GEO_OPT"
            print "    TYPE MINIMIZATION"
            print "    KEEP_SPACE_GROUP F"
            print "    OPTIMIZER BFGS"
            print "    &BFGS"
            print "      TRUST_RADIUS 0.2"
            print "    &END BFGS"
            print "    MAX_ITER 500"
            print "    MAX_DR 3E-3"
            print "    RMS_DR 1.5E-3"
            print "    MAX_FORCE 4.5E-4"
            print "    RMS_FORCE 3E-4"
            print "  &END GEO_OPT"
            print "  &PRINT"
            print "    &TRAJECTORY"
            print "      FORMAT xyz"
            print "    &END TRAJECTORY"
            print "  &END PRINT"
            print "&END MOTION"
        }
        ' "$CIF2CELL_IN" > temp && mv temp "$CIF2CELL_IN"

    # Using xTB method
    read -p "Using xTB method ? [Y/n] " cp2k_xtb; cp2k_xtb=${cp2k_xtb:-"y"}

    if [[ "$cp2k_xtb" == [Yy] ]]; then
        sed -i -e '/^\s*&KIND/,/^\s*&END KIND/d' \
            -e '/^\s*&XC$/,/^\s*&END MGRID$/d' \
            -e '/^\s*&DIAGONALIZATION$/,/^\s*&END MIXING$/d' \
            -e '/BASIS_SET_FILE_NAME/d' \
            -e '/POTENTIAL_FILE_NAME/d' \
            "$CIF2CELL_IN"

        sed -i '/^\s*&QS$/,/^\s*&END QS$/ {
            /^\s*&END QS$/i\
        METHOD xTB\
        &xTB\
          DO_EWALD T\
          CHECK_ATOMIC_CHARGES F\
          &PARAMETER\
          DISPERSION_PARAMETER_FILE dftd3.dat\
          PARAM_FILE_NAME xTB_parameters\
          &END PARAMETER\
        &END xTB
            }' "$CIF2CELL_IN"

        sed -i '/^\s*&SCF$/,/^\s*&END SCF$/ {
            /^\s*&END SCF$/i\
        &OT\
          PRECONDITIONER FULL_SINGLE_INVERSE\
          MINIMIZER DIIS\
          LINESEARCH 2PNT\
          ALGORITHM STRICT\
        &END OT\
        &OUTER_SCF\
          MAX_SCF 20\
          EPS_SCF 1.0E-06\
        &END OUTER_SCF
            }' "$CIF2CELL_IN"
    fi

    # Using vc-relax
    read -p "Relax the cell ? [y/N] " cp2k_cell_opt; cp2k_cell_opt=${cp2k_cell_opt:-"n"} 

    if [[ "$cp2k_cell_opt" == [Yy] ]]; then
        cp2k_CALC_TYPE="CELL_OPT"
        sed -i -e 's/GEO_OPT/CELL_OPT/g' \
            -e 's/FORMAT xyz/FORMAT pdb/g' \
            -e 's/MINIMIZATION/DIRECT_CELL_OPT/g' \
            "$CIF2CELL_IN"

        sed -i '/^\s*&CELL_OPT$/,/^\s*&END CELL_OPT$/ {
            /^\s*&END CELL_OPT$/i\
    EXTERNAL_PRESSURE  1.01325E+00\
    PRESSURE_TOLERANCE 100\
    CONSTRAINT NONE\
    KEEP_ANGLES F\
    KEEP_SYMMETRY F
            }' "$CIF2CELL_IN"

        sed -i '/^\s*&FORCE_EVAL$/,/^\s*&END FORCE_EVAL$/ {
            /^\s*&END FORCE_EVAL$/i\
  &PRINT\
    &STRESS_TENSOR ON\
    &END STRESS_TENSOR\
  &END PRINT\
  STRESS_TENSOR ANALYTICAL
            }' "$CIF2CELL_IN"

        sed -i '/^ * &END TRAJECTORY$/a\
    &RESTART\
      BACKUP_COPIES 0\
    &END RESTART\
    &RESTART_HISTORY OFF\
    &END RESTART_HISTORY' "$CIF2CELL_IN"
    fi
    
    mv $CIF2CELL_IN ${CIF_FILE%.*}.$cp2k_CALC_TYPE.inp

}

in2cif() {
    
    IN_FILE="$1"

    awk '
        /A *=/ { A = $NF }
        /CELL_PARAMETERS *{alat}/ {
            getline
            print "%BLOCK LATTICE_CART"
            for (i = 0; i < 3; i++) {
                split($0, a)
                for (j = 1; j <= 3; j++)
                    printf "%.15f%s", a[j]*A, (j==3 ? "\n" : "  ")
                getline
            }
            print "%ENDBLOCK LATTICE_CART"
            print ""
        }

        /ATOMIC_POSITIONS *{crystal}/ {
            print "%BLOCK POSITIONS_FRAC"
            while (getline) {
                if ($2 !~ /^[0-9.+-]/) break
                print
            }
            print "%ENDBLOCK POSITIONS_FRAC"
            exit
        }
    ' "$IN_FILE" > $IN_FILE.cell
        
    PYTHONWARNINGS="ignore" ase convert $IN_FILE.cell $IN_FILE.cif
}

in2conv () {

    IN_FILE="$1"
    CONVTEST=${CONVTEST:-"ecut"}

    read -p "Set test values for $CONVTEST test: [40 50 60 70 80] " TEST_VALUES ; TEST_VALUES=${TEST_VALUES:-"40 50 60 70 80"}
    mkdir -p ./${IN_FILE%.*}.$CONVTEST

    for TEST_VALUE in ${TEST_VALUES} ; do
        IN_FILE_OUT="./${IN_FILE%.*}.$CONVTEST/${IN_FILE%.*}.conv_$TEST_VALUE.in"
        cp $IN_FILE $IN_FILE_OUT
        if grep -q "^ *ecutwfc *= *" "$IN_FILE_OUT"; then
            sed -i "s/^ *ecutwfc *= .*/  ecutwfc = $TEST_VALUE/" "$IN_FILE_OUT"
        else
            sed -i "/ntyp/a\  ecutwfc = $TEST_VALUE" "$IN_FILE_OUT"
        fi 
    done
}

out_scf () {

    OUT_FILE="$1"
    
    echo "Cycle,Energy,Accuracy" > "${OUT_FILE%.*}_scf.csv"
    paste -d ' ' \
        <(grep  "iteration #" "$OUT_FILE" | awk -F'[#]' '{printf "%d\n", $2}') \
        <(grep  "total energy              =" "$OUT_FILE" | awk -F'[=]' '{printf "%.10f\n", $2}') \
        <(grep  "estimated scf accuracy" "$OUT_FILE" | awk -F'[<]' '{printf "%.6e\n", $2}') \
    | awk '{printf "%d,%.10f,%.6e\n", $1, $2, $3}' >> "${OUT_FILE%.*}_scf.csv"

}

out_relax () {

    OUT_FILE="$1"

    paste -d ' ' \
        <(grep -e "!    total energy" "$OUT_FILE" | awk -F'=' '{printf "%.10f\n", $2}') \
        <(grep "Energy error" "$OUT_FILE" | awk '{print $(NF-1)}') \
        <(grep "Gradient error" "$OUT_FILE" | awk '{print $(NF-1)}') \
        <(grep "Cell gradient error" "$OUT_FILE" | awk '{print $(NF-1)}') \
    | awk -v cell_relax="$(grep -c "Cell gradient error" "$OUT_FILE")" \
        -v scf_num="$(grep -c "Energy error" "$OUT_FILE")" \
        'BEGIN {
            if (cell_relax == 0) {
                print "Cycle,Energy,deltaE,deltaGrad"
            } else {
                print "Cycle,Energy,deltaE,deltaGrad,deltaCell"
            }
        }
        NR <= scf_num {
            if (cell_relax == 0) {
                printf "%d,%.10f,%s,%s\n", NR, $1, $2, $3
            } else {
                printf "%d,%.10f,%s,%s,%s\n", NR, $1, $2, $3, $4
            }
        }' > "${OUT_FILE%.*}_grad.csv"
        
    ase convert -i espresso-out -o cif -n 0 "$OUT_FILE" "${OUT_FILE%.*}_start.cif"
    ase convert -i espresso-out -o cif -n -1 "$OUT_FILE" "${OUT_FILE%.*}_final.cif"
    ase convert -i espresso-out -o extxyz "$OUT_FILE" "${OUT_FILE%.*}_trj.xyz"
}

cif_gen () {
    echo "Choose lattice type:"
    echo "1. Cubic   2. Tetragonal   3. Orthorhombic"
    echo "4. Hexagonal   5. Monoclinic   6. Triclinic"
    read -p "Choice (1-6): " cif_gen_lattice

    # Default cell
    a=20.0; b=20.0; c=20.0
    alpha=90.0; beta=90.0; gamma=90.0

    case $cif_gen_lattice in
        1) read -p "a = " a; b=$a; c=$a ;;
        2) read -p "a = " a; read -p "c = " c; b=$a ;;
        3) read -p "a = " a; read -p "b = " b; read -p "c = " c ;;
        4) read -p "a = " a; b=$a; read -p "c = " c; gamma=120.0 ;;
        5) read -p "a = " a; read -p "b = " b; read -p "c = " c; read -p "β = " beta ;;
        6) read -p "a = " a; read -p "b = " b; read -p "c = " c
           read -p "alpha = " alpha; read -p "beta = " beta; read -p "gamma = " gamma ;;
    esac

    read -p "Space group (e.g. P1): " space_group; space_group=${space_group:-"P1"}

    echo "Enter atoms (Symbol x y z), type 'done' to finish:"
    atom_lines=()
    while true; do
        read -p "> " line
        [[ "$line" == "done" ]] && break
        atom_lines+=("$line")
    done

    CIF_GEN_FILE="${1:-structure}.cif"

    awk -v a="$a" -v b="$b" -v c="$c" \
        -v alpha="$alpha" -v beta="$beta" -v gamma="$gamma" \
        -v spg="$space_group" '
        BEGIN {
            print "data_cifgen"
            printf "_cell_length_a            %s\n", a
            printf "_cell_length_b            %s\n", b
            printf "_cell_length_c            %s\n", c
            printf "_cell_angle_alpha         %s\n", alpha
            printf "_cell_angle_beta          %s\n", beta
            printf "_cell_angle_gamma         %s\n", gamma
            printf "_space_group_name_H-M_alt '\''%s'\''\n", spg
            print "\nloop_"
            print "  _atom_site_type_symbol"
            print "  _atom_site_fract_x"
            print "  _atom_site_fract_y"
            print "  _atom_site_fract_z"
            print "  _atom_site_occupancy"
        }
        {
            printf "  %-4s %12s %12s %12s     1.000\n", $1, $2, $3, $4
        }
    ' <<< "$(printf "%s\n" "${atom_lines[@]}")" > "$CIF_GEN_FILE"

    echo "CIF file saved: $CIF_GEN_FILE"
}

############################################# main ##########################################
READ_FILE="$1"

# Validate file input or prompt if invalid
if [[ ! -f "$READ_FILE" || ! -r "$READ_FILE" ]]; then
    echo "Input file '$READ_FILE' is missing or unreadable."
    read -p "Do you want to generate a CIF file muanally ? [Y/n]: " cif_gen_choice; cif_gen_choice=${cif_gen_choice:-Y}
    
    if [[ "$cif_gen_choice" =~ ^[Yy]$ ]]; then
        cif_gen $READ_FILE
        exit 0
    else
        echo "Exiting script. Please check the filename"
        exit 1
    fi
fi

echo "Successfully loaded input file: $READ_FILE"

# Detect file type
EXT="${READ_FILE##*.}"
case "$EXT" in
    cif) 
        echo "Detected file type: CIF"
        CIF_FILE=$READ_FILE
        APP_FORMAT="${2:-qe}"
        case "$APP_FORMAT" in
            qe) 
                CALC_TYPE="${3:-scf}"
                echo "Starting to prepare the $CALC_TYPE input files for Quantum ESPRESSO program"
                cif2qe $CIF_FILE $CALC_TYPE
                ;;
            vasp) 
                echo "Starting to prepare the input files for VASP program"
                cif2vasp $CIF_FILE
                ;;
            cp2k) 
                echo "Starting to prepare the input files for CP2K program"
                cif2cp2k $CIF_FILE
                ;;
            *) 
                echo "Unknown application format: $APP_FORMAT"
                exit 1
                ;;
        esac
        ;;
    in)  
        echo "Detected file type: Quantum ESPRESSO IN"; 
        IN_FILE=$READ_FILE 
        EXP_FORMAT="${2:-cif}"
        case "$EXP_FORMAT" in
            cif) 
                echo "Converting the Quantum ESPRESSO input files to CIF"
                in2cif $IN_FILE
                ;;
            conv) 
                echo "Starting to prepare conv test files for Quantum ESPRESSO program"
                in2conv $IN_FILE
                ;;
            *) 
                echo "Unknown format: $EXP_FORMAT"
                exit 1
                ;;
        esac
        ;;
    out) 
        echo "Detected file type: Quantum ESPRESSO OUT"; 
        OUT_FILE=$READ_FILE 
        OUT_ANALYSIS="${2:-relax}"
        case "$OUT_ANALYSIS" in
            scf) 
                echo "Analysis the Quantum ESPRESSO scf output files"
                out_scf $OUT_FILE
                ;;
            relax) 
                echo "Analysis the Quantum ESPRESSO relax or vc-relax output files"
                out_relax $OUT_FILE
                ;;
            *) 
                echo "Unknown analysis type: $OUT_ANALYSIS"
                exit 1
                ;;
        esac
        ;;
    *)   
        echo "Unknown file extension: .$EXT" 
        ;;
esac
