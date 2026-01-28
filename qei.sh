#!/bin/bash
# ==============================================================================
# Script Name: qei.sh
# Description: Quantum ESPRESSO Input Generation & Conversion Tool
#
# Summary:
#   A multi-functional pre-processing and analysis tool. It handles:
#   1. Structure -> Input: Converts CIF files to inputs for QE, VASP, or CP2K.
#   2. Input -> Structure: Converts QE input files back to CIF.
#   3. Input -> Conv Test: Generates convergence test files from a QE input.
#   4. Output -> Analysis: Extracts key data from QE output files.
#
# Usage:
#   qei.sh <file.cif> [format] [type]
#   qei.sh <file.in>  [action]
#   qei.sh <file.out> [type]
#
# Version: 0.4.1 (2026-02-28)
# ==============================================================================

# ==============================================================================
# Configuration & Path Definitions
# ==============================================================================

SCRIPT_NAME=$(basename "$0")
VERSION="0.4.1"
DATE="2026-02-28"

ASE_PATH="ase"
CIF2CELL_PATH="cif2cell"
QE_PSEUDOLIB="/opt/_pseudo/upf_PBEsol"
VASP_PSEUDOLIB="/opt/_pseudo/paw64_PBE"
CP2K_DATA_DIR=/opt/cp2k-2025.2/data

TMP_DIR="./tmp"
PSEUDO_DIR="./pseudo"
CIF2CELL_IN="cif2cell.in"

# ---------- Help Function ----------
print_help() {
  cat <<EOF
${SCRIPT_NAME} — Input Generation & Conversion v${VERSION} (${DATE})

USAGE
  1. From Structure (CIF):
     ${SCRIPT_NAME} <file.cif> [qe|vasp|cp2k] [calc_type]
       • Generates input files for the specified code.
       • Default format: qe
       • Default calc_type: scf
       • Example: ${SCRIPT_NAME} structure.cif qe relax

  2. From Input (.in):
     ${SCRIPT_NAME} <file.in> [cif|conv]
       • cif  : Converts QE input back to CIF format (default).
       • conv : Generates convergence test files (k-points/cutoff).
       • Example: ${SCRIPT_NAME} scf.in conv

  3. From Output (.out):
     ${SCRIPT_NAME} <file.out> [scf|relax|vc-relax]
       • Analyzes output logs.
       • Default analysis: relax
       • Example: ${SCRIPT_NAME} relax.out vc-relax

  4. Interactive / Manual:
     ${SCRIPT_NAME} <missing_file>
       • Triggers manual CIF generation wizard.

FLAGS
  -h, --help    Show this help message.

EOF
}

# parse z_valence from UPF
_parse_zval() {
  local upf="$1"; local z
  if [[ "$upf" == *-gbrv* ]]; then
    z=$(grep -m1 "Z valence" "$upf" | awk '{printf "%.1f", $1}')
  else
    z=$(grep -oP -m1 'z_valence\s*=\s*["'\'']?\s*\K[0-9.]+' "$upf" | awk '{printf "%.1f", $1}')
  fi
  [[ -n "$z" ]] || z=0
  printf '%s\n' "$z"
}

# count atoms per element by parsing ATOMIC_POSITIONS block
_count_atoms() {
  local infile="$1"; shift
  awk '
    /^ATOMIC_POSITIONS/ {inpos=1; next}
    inpos && NF==0 {inpos=0}
    inpos {count[$1]++}
    END{for(e in count) printf "%s:%d\n", e, count[e]}
  ' "$infile" | sort
}

# awk in-place via temp, POSIX-safe
_awk() {
  local tmp file
  local args=("$@")
  file="${args[-1]}"    
  unset "args[-1]"      
  tmp=$(mktemp "${file##*/}.XXXX") || return 1
  awk "${args[@]}" "$file" >"$tmp" && mv "$tmp" "$file"
}

cif2qe() {

    CIF_FILE="$1"
    CALC_TYPE="$2"
    
    PP_FILE_NAME="-pseudo.UPF" 
    SCF_K_RESOLUTION="0.35" 
    DEGAUSS="0.001"
    $CIF2CELL_PATH "$CIF_FILE" -p pwscf -o $CIF2CELL_IN \
        --setup-all --no-reduce --print-digits=10 --k-resolution=$SCF_K_RESOLUTION --pwscf-pseudostring=$PP_FILE_NAME 

    # Prepare ELEM_INFO and PSEUDOs
    read -p "==> Assigin pseudo type, DFT functional, and ecutwfc value [US/paw/nc PBEsol 50 Ry]: " PSEUDO_TYPE DFTINPUT ECUT; PSEUDO_TYPE=${PSEUDO_TYPE:-us}; DFTINPUT=${DFTINPUT:-pbesol}; ECUT=${ECUT:-50}
    read -p "==> Set charge and magnetization [0 0]: " CHARGE MAG; CHARGE=${CHARGE:-0}; MAG=${MAG:-n}
    total_electrons=$(echo "0.0 - $CHARGE" | bc -l)
    declare -A ELEM_INFO
    declare -A U_PRESETS=(
        ["Ti"]="3d 4.0" 
        ["Fe"]="3d 5.0" 
        ["Cu"]="3d 4.0" 
        ["Ni"]="3d 6.0" 
        ["Co"]="3d 5.0" 
        ["V"]="3d 3.0"  
        ["Cr"]="3d 3.0" 
        ["Mn"]="3d 4.0" 
        ["Zn"]="3d 8.0" 
        ["Sc"]="3d 2.0" 
    )

    while IFS=':' read -r elem count; do
        ELEM_INFO["$elem"]="$count"
    done < <(_count_atoms "$CIF2CELL_IN")

    read -p "==> Copy pseudo files to ${PSEUDO_DIR} ? [y/N] " CP_PSEUDO; CP_PSEUDO=${CP_PSEUDO:-n}
    for ELEM in "${!ELEM_INFO[@]}"; do
        files=($(ls $QE_PSEUDOLIB/${ELEM}-*${DFTINPUT}-*${PSEUDO_TYPE}*.UPF 2>/dev/null))
        if [ ${#files[@]} -eq 0 ]; then
            echo "No pseudopotential files matching current criteria"
            exit 1
        elif [ ${#files[@]} -eq 1 ]; then
            PSEUDO_NEW=$(basename "${files[0]}")
        else
            echo "---------------------------------------------------------------"
            echo "Setup pseudo file for $ELEM :"
            echo "---------------------------------------------------------------"
            for i in "${!files[@]}"; do
                echo "    $i: $(basename "${files[$i]}")"
            done
            read -p "==> Setup $ELEM with [0] " PSEUDO_CHOICE
            PSEUDO_CHOICE=$([[ "$PSEUDO_CHOICE" =~ ^[0-9]+$ ]] && ((PSEUDO_CHOICE >= 0 && PSEUDO_CHOICE < ${#files[@]})) && echo "$PSEUDO_CHOICE" || echo 0)
            PSEUDO_NEW=$(basename "${files[$PSEUDO_CHOICE]}")
        fi

        if [[ "$CP_PSEUDO" =~ ^[Yy]$ ]]; then
            mkdir -p "$PSEUDO_DIR"
            cp "$QE_PSEUDOLIB/$PSEUDO_NEW" "$PSEUDO_DIR"
        fi

        z_val=$(_parse_zval "$QE_PSEUDOLIB/$PSEUDO_NEW")
        is_ts=$([[ -v U_PRESETS[$ELEM] ]] && echo true || echo false)
        total_electrons=$(echo "$total_electrons + $z_val * ${ELEM_INFO["$ELEM"]}" | bc -l)
        ELEM_INFO["$ELEM"]="$PSEUDO_NEW|$z_val|$is_ts|${ELEM_INFO["$ELEM"]}"
        sed -i "s/$(grep -o '\S*\.UPF' "$CIF2CELL_IN" | grep "${ELEM}-")/${PSEUDO_NEW}/g" "$CIF2CELL_IN"
    done

    # Verify final ELEM_INFO
    hoco_num=$(echo "$total_electrons / 2" | bc -l | xargs printf "%.0f")
    luco_num=$((hoco_num + 1))
    co_min=$((hoco_num - 2))
    co_max=$((hoco_num + 3))
    suggested_band=$(awk -v ne="$total_electrons" 'BEGIN { printf "%d", int(ne/2*1.25) + 4 }')
    echo -e "\nFinal PSEUDO info:"
    echo "---------------------------------------------------------------"
    printf "%-5s %-25s %-8s %-10s %s\n" "Elem" "Pseudo File" "z_val" "is_ts" "Atom Count"
    echo "---------------------------------------------------------------"
    for elem in "${!ELEM_INFO[@]}"; do
        IFS='|' read -r pseudo z_val is_ts count <<< "${ELEM_INFO[$elem]}"
        printf "%-5s %-25s %-8s %-10s %d\n" "$elem" "$pseudo" "$z_val" "$is_ts" "$count"
    done
    echo "---------------------------------------------------------------"
    echo "Nele: $(printf "%.1f" "$total_electrons"); Nbnds: $suggested_band ($co_min → $co_max); HOCO: $hoco_num; LUCO: $luco_num"





    # Set up Huaabrd U
    current_ts_metals=()  
    for ELEM in "${!ELEM_INFO[@]}"; do
        if [ "$(echo "${ELEM_INFO[$ELEM]}" | cut -d'|' -f3)" = "true" ]; then
            current_ts_metals+=("$ELEM")  
        fi
    done

    if [ ${#current_ts_metals[@]} -gt 0 ]; then
        read -p "==> Using Hubbard U for transition metal: ${current_ts_metals[*]}? [y/N] " SET_HUBBARDU; SET_HUBBARDU=${SET_HUBBARDU:-N}
        if [[ "$SET_HUBBARDU" =~ ^[Yy]$ ]]; then
            echo "" >> "$CIF2CELL_IN"
            echo "HUBBARD {ortho-atomic}" >> "$CIF2CELL_IN"

            for ELEM in "${current_ts_metals[@]}"; do
                if [ -n "${U_PRESETS[$ELEM]}" ]; then
                    orbital=${U_PRESETS[$ELEM]% *}  # Extract orbital (e.g., 3d)
                    u_value=${U_PRESETS[$ELEM]#* }  # Extract U value (e.g., 4.0)
                    echo "U     $ELEM-$orbital  $u_value" >> "$CIF2CELL_IN"
                    echo "Automatically set U value for $ELEM-$orbital to preset: $u_value eV"
                else
                    # Ask for input if no preset exists
                    read -p "No preset U value found for $ELEM, please enter orbital (e.g., 3d/4d): " orbital
                    read -p "Please enter U value for $ELEM-$orbital (eV): " u_value
                    # Validate input is a number
                    while ! [[ "$u_value" =~ ^[0-9]+(\.[0-9]+)?$ ]]; do
                        echo "Input error! U value must be a number (e.g., 4.0)"
                        read -p "Please re-enter U value for $ELEM-$orbital (eV): " u_value
                    done
                    echo "U     $ELEM-$orbital  $u_value" >> "$CIF2CELL_IN"
                    echo "Set U value for $ELEM-$orbital to: $u_value eV"
                fi
            done
        fi
    fi





    # &CONTROL, &SYSTEM, &ELECTRONS, &IONS, &CELL and K_POINTS blocks
    _awk -v calc_type="$CALC_TYPE" -v charge="$CHARGE" -v magnetization="$MAG" -v cp_pseudo="$CP_PSEUDO"\
        -v pseudo_dir="$PSEUDO_DIR" -v prefix="${CIF_FILE%.*}" -v tmp_dir="$TMP_DIR" -v num_band="$suggested_band"\
        -v dftinput="$DFTINPUT" -v ecut="$ECUT" -v pseudo_type="$PSEUDO_TYPE" \
        'BEGIN {
            in_system = 0; done = 0 
            print "&CONTROL"
            print "  calculation = '\''" calc_type "'\''"
            print "  prefix = '\''" prefix "'\''"
            if (calc_type == "scf") {
                print "  verbosity = '\''high'\''"
            }
            print "  outdir = '\''" tmp_dir "'\''"
            if ( cp_pseudo ~ /^[Yy]$/) {
                print "  pseudo_dir = '\''" pseudo_dir "'\''"
            }
            print "  !restart_mode = '\''restart'\''"
            if (calc_type == "relax" || calc_type == "vc-relax") {
                print "  nstep = 200"
                print "  forc_conv_thr = 1.0d-3"
                print "  etot_conv_thr = 1.0d-4"
            }
            print "/"
        }
        /^&SYSTEM/ { in_system = 1 }
        in_system && /^\/$/ && !done {
            if (calc_type != "relax" || calc_type != "vc-relax") {
                print "  nbnd = " num_band
            }
            print "  ecutwfc = " ecut
            if (pseudo_type == "us") {
                print "  ecutrho = " (ecut * 10)
            } else if (pseudo_type == "paw") {
                print "  ecutrho = " (ecut * 6)
            } else {
                print "  ecutrho = " (ecut * 4)
            }
            
            if (charge != "0") print "  tot_charge = " charge
            if (magnetization != "n") {
                print "  nspin = 2"
                print "  tot_magnetization = " magnetization
                print "  !starting_magnetization(1) = 0.8"
                print "  occupations = '\''smearing'\''"
                print "  smearing = '\''mv'\''"
                print "  degauss = 0.001"
            }

            if (dftinput != "PBEsol") {
                print "  input_dft = '\''" (dftinput == "hsesol" ? "sla+pw+hse+psc" : dftinput) "'\''"
            }
            print "  vdw_corr = '\''dft-d3'\''"
            print "  dftd3_version = 4"
            done = 1
        }
        /^\/$/ {
            print
            if (in_system) {
                print "&ELECTRONS"
                print "  mixing_beta = 0.7"
                print "  conv_thr = 1.0d-6"
                print "  electron_maxstep = 200"
                print "  !diagonalization = '\''cg'\''"
                print "  !startingwfc = '\''random'\''"
                if (calc_type == "relax" || calc_type == "vc-relax") {
                    print "  scf_must_converge = .false."
                }
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
        NR > 9 { print }
        ' "$CIF2CELL_IN"

    mv $CIF2CELL_IN ${CIF_FILE%.*}.$CALC_TYPE.in
}  

cif2vasp () {
    
    CIF_FILE="$1"
    CALC_TYPE="$2"
    
    mkdir -p ${CIF_FILE%.*}.$CALC_TYPE && cd ${CIF_FILE%.*}.$CALC_TYPE

    VASP_PAWLIB="$VASP_PSEUDOLIB" VASP_PP_PRIORITY="_d,_pv,_sv,_h,_s" $CIF2CELL_PATH "../$CIF_FILE" -p vasp --vasp-format=5\
        --setup-all --no-reduce --print-digits=10 --k-resolution=$SCF_K_RESOLUTION  
    NBANDS=$(grep NBANDS INCAR | awk '{print $NF}')
    
    
    direct_line=$(grep -n "^Direct" POSCAR | cut -d: -f1)
    types_line=$((direct_line - 2))  # Line with element symbols (e.g., "Ti   O")
    counts_line=$((direct_line - 1)) # Line with atom counts (e.g., "2   4")
    atom_types=($(sed -n "${types_line}p" POSCAR | tr -s '[:space:]' '\n' | grep -v '^$'))
    atom_counts=($(sed -n "${counts_line}p" POSCAR | tr -s '[:space:]' '\n' | grep -v '^$' | awk '{print int($1)}'))

    declare -A valence_electrons
    # Store data in an associative array (key: element, value: count)
    for i in "${!atom_types[@]}"; do
        element="${atom_types[$i]}"
        
        # Find the "PAW_" line in POTCAR for this element (with leading space to avoid partial matches)
        # Example: Matches " PAW_PBE Ti_pv ..." for element "Ti"
        paw_line=$(grep -n "^ *PAW_.* $element" POTCAR | cut -d: -f1)
        
        if [ -z "$paw_line" ]; then
            echo "Error: No POTCAR entry found for element '$element' (searched for 'PAW_... $element')"
            exit 1
        fi
        
        # Valence electrons are in the line immediately after the "PAW_" line
        valence_line=$((paw_line + 1))
        valence=$(sed -n "${valence_line}p" POTCAR | awk '{print int($1)}')  # Convert to integer
        
        if [ -z "$valence" ] || ! [[ "$valence" =~ ^[0-9]+$ ]]; then
            echo "Error: Failed to extract valence electrons for '$element' from POTCAR line $valence_line"
            exit 1
        fi
        
        valence_electrons["$element"]=$valence
    done

    total_nelect=0

    # Print header
    printf "%-8s %6s %12s %12s\n" "Element" "Count" "Valence e⁻" "Total e⁻"

    # Print single separator line (45 dashes to cover the header width)
    printf "%0.s-" {1..40}; echo

    # Print element data
    for i in "${!atom_types[@]}"; do
        element="${atom_types[$i]}"
        count="${atom_counts[$i]}"
        valence="${valence_electrons[$element]}"
        per_element_total=$((count * valence))
        total_nelect=$((total_nelect + per_element_total))
        
        # Align columns with header
        printf "%-8s %6d %10d %10d\n" "$element" "$count" "$valence" "$per_element_total"
    done

    # Print second separator line
    printf "%0.s-" {1..40}; echo

    # Print total
    printf "Total electrons (NELECT) = %d\n" "$total_nelect"
    read -p "==> Set charge and magnetization [0 0]: " -a input
    # Extract CHARGE (first element) with default 0
    CHARGE=${input[0]:-0}

    # Extract MAG (remaining elements) with default "n", joined by spaces
    MAG="${input[@]:1}"  # All elements from index 1 onwards
    MAG=${MAG:-n}        # Default to "n" if no input for MAG

    NUM_ELECT=$(echo "$total_nelect - $CHARGE" | bc)
    awk -v nbands="$NBANDS" -v calc="$CALC_TYPE" -v nelect="$NUM_ELECT" \
        -v charge="$CHARGE" -v magnetization="$MAG" \
        'BEGIN {
            print "&CONTROL"
            print "# ISTART = 0 (start from scratch)"
            print "LREAL = Auto"
            print "NCORE = 4 (number of cores per band)"
            print "KPAR = 3 (Divides k-grid into separate groups)"
            print ""
            print "&ELECTRONS"
            print "GGA = PS (PBEsol functional)"
            print "ENCUT = 450 (plane-wave cutoff, eV)"
            print "IVDW = 12 (D3 van der Waals correction)"
            if (charge != "0") print "NELECT = " nelect
            if (magnetization != "n") {
                print "ISPIN = 2"
                print "MAGMOM = " magnetization
            }
            print "NELM = 200 (max electronic steps)"
            print "NELMIN = 6 (min electronic steps)"
            print "EDIFF = 1E-08 (electronic convergence, eV)"
            print "ISMEAR = -5 (Gaussian smearing for DOS)"
            print "SIGMA = 0.05 (smearing width, eV)"
            print ""

        # Set ionic relaxation parameters based on calculation type
        if (calc == "scf") {
            print "IBRION = -1 (no ionic relaxation for SCF)"
            print "NSW = 0 (no ionic steps for SCF)"
            print "PREC = Accurate"
            print "# ADDGRID = .TRUE."
            print "# LASPH = .TRUE."
            print "# ICHARG = 11"
            print "NBANDS = " nbands
            print "NEDOS = 2001 (DOS energy grid points)"
            print "LORBIT = 11 (PAW radii for projected DOS)"
            print "# LVHAR = .TRUE. (total potential output)"
        } else {
            print "&IONS"
            if (calc == "relax") {print "ISIF = 2"}
            if (calc == "vc-relax") {print "ISIF = 3"}
            print "IBRION = 1 (0-MD, 1-Quasi-New, 2-CG)"
            print "POTIM = 0.5 (step width in ionic relaxations)"
            print "NSW = 200 (max ionic steps)"
            print "EDIFFG = -0.01"
            print "NWRITE = 1 (low-level output)"
        }
        }' > INCAR && cd ..
}

cif2cp2k () {

    CIF_FILE="$1"
    CALC_TYPE="$2"
    $CIF2CELL_PATH "$CIF_FILE" -p cp2k -o $CIF2CELL_IN --no-reduce --print-digits=10

    case "$CALC_TYPE" in
    scf)        CALC_TYPE="ENERGY" ;;
    relax)      CALC_TYPE="GEO_OPT" ;;
    vc-relax)   CALC_TYPE="CELL_OPT" ;;
    freq)       CALC_TYPE="VIBRATIONAL_ANALYSIS" ;;
    esac

    read -p "==> Set charge and magnetization [0 0]: " CHARGE MAG; CHARGE=${CHARGE:-0}; MAG=${MAG:-n}
    read -p "==> Using xTB method ? [Y/n] " cp2k_xtb; cp2k_xtb=${cp2k_xtb:-"y"}
    _awk -v proj="${CIF_FILE%.*}" -v calc_type="$CALC_TYPE" -v xtb_input="$cp2k_xtb" -v charge="$CHARGE" -v magnetization="$MAG"\
        'BEGIN {
            print "&GLOBAL"
            print "  PROJECT " proj
            if (calc_type == "GEO_OPT" || calc_type == "CELL_OPT") {
                print "  PRINT_LEVEL LOW"
            }
            print "  RUN_TYPE " calc_type
            print "&END GLOBAL"
            print ""
            print "&FORCE_EVAL"
            print "  METHOD Quickstep"
            print "  &SUBSYS"
            in_coord = 0
        }

        /&COORD/     { in_coord = 1 }
        /&END COORD/ { in_coord = 0 }

        NR > 9 {
            if (in_coord && NF == 4 && $1 ~ /^[A-Za-z]+$/) {
                elem = (length($1) == 1) ? $1 " " : $1
                printf("      %-2s %20.15f %20.15f %20.15f\n", elem, $2, $3, $4)
                if (!(elems_seen[$1]++)) elems[++n] = $1
            } else {
                if (NF > 0) print "    " $0
            }
        }

        END {
            if (xtb_input !~ /^[Yy]$/) {
                for (i = 1; i <= n; i++) {
                    printf "    &KIND %s\n", elems[i]
                    printf "      ELEMENT %s\n", elems[i]
                    print  "      BASIS_SET DZVP-MOLOPT-SR-GTH"
                    print  "      POTENTIAL GTH-PBE"
                    if (magnetization != "n") {
                        print  "      # MAGNETIZATION 0.5"
                    }
                    print  "    &END KIND"
                }
            }
            print "  &END SUBSYS"
            print ""
            print "  &DFT"
            if (xtb_input !~ /^[Yy]$/) {
                print "    BASIS_SET_FILE_NAME  BASIS_MOLOPT"
                print "    POTENTIAL_FILE_NAME  GTH_POTENTIALS"
            }
            print "    CHARGE    " charge
            if (magnetization != "n") {
                print "    MULTIPLICITY    " (magnetization + 1)
                print "    UKS TRUE"
            }
            if (xtb_input !~ /^[Yy]$/) {
                print "    &KPOINTS"
                print "      SCHEME GAMMA"
                print "    &END KPOINTS"
            }
            print "    &QS"
            print "      EPS_DEFAULT 1.0E-12"
            if (xtb_input ~ /^[Yy]$/) {
                print "      METHOD xTB"
                print "      &xTB"
                print "        DO_EWALD T"
                print "        CHECK_ATOMIC_CHARGES F"
                print "      &END xTB"
            }
            print "    &END QS"
            if (xtb_input !~ /^[Yy]$/) {
                print "    &XC"
                print "      &XC_FUNCTIONAL"
                print "        &PBE"
                print "          PARAMETRIZATION PBESOL"
                print "        &END PBE"
                print "      &END XC_FUNCTIONAL"
                print "      &VDW_POTENTIAL"
                print "        POTENTIAL_TYPE PAIR_POTENTIAL"
                print "        &PAIR_POTENTIAL"
                print "          TYPE DFTD3(BJ)"
                print "          REFERENCE_FUNCTIONAL PBEsol"
                print "        &END PAIR_POTENTIAL"
                print "      &END VDW_POTENTIAL"
                print "    &END XC"
                print "    &MGRID"
                print "      CUTOFF  400"
                print "      REL_CUTOFF  60"
                print "    &END MGRID"
            }
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
            if (xtb_input !~ /^[Yy]$/) {
                print "      &DIAGONALIZATION"
                print "        ALGORITHM STANDARD"
                print "      &END DIAGONALIZATION"
                print "      &MIXING"
                print "        METHOD BROYDEN_MIXING"
                print "        ALPHA 0.4"
                print "        NBROYDEN 8"
                print "      &END MIXING"
            } else {
                print "      &OT"
                print "        PRECONDITIONER FULL_SINGLE_INVERSE"
                print "        MINIMIZER DIIS"
                print "        LINESEARCH 2PNT"
                print "        ALGORITHM STRICT"
                print "      &END OT"
                print "      &OUTER_SCF"
                print "        MAX_SCF 20"
                print "        EPS_SCF 1.0E-06"
                print "      &END OUTER_SCF"
            }
            print "    &END SCF"
            print "  &END DFT"
            if (calc_type == "CELL_OPT") {
                print "  &PRINT"
                print "    &STRESS_TENSOR ON"
                print "    &END STRESS_TENSOR"
                print "  &END PRINT"
                print "  STRESS_TENSOR ANALYTICAL"
            }
            print "&END FORCE_EVAL"
            if (calc_type == "GEO_OPT" || calc_type == "CELL_OPT") {
                print ""
                print "&MOTION"
                if (calc_type == "GEO_OPT") {
                    print "  &GEO_OPT"
                    print "    TYPE MINIMIZATION"
                } else {
                    print "  &CELL_OPT"
                    print "    EXTERNAL_PRESSURE  1.01325E+00"
                    print "    PRESSURE_TOLERANCE 100"
                    print "    CONSTRAINT NONE"
                    print "    KEEP_ANGLES F"
                    print "    KEEP_SYMMETRY F"
                }
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
                if (calc_type == "GEO_OPT") {
                    print "  &END GEO_OPT"
                } else {
                    print "  &END CELL_OPT"
                }
                print "  &PRINT"
                print "    &TRAJECTORY"
                print "      FORMAT pdb"
                print "    &END TRAJECTORY"
                print "    &RESTART"
                print "      BACKUP_COPIES 0"
                print "    &END RESTART"
                print "    &RESTART_HISTORY OFF"
                print "    &END RESTART_HISTORY"
                print "  &END PRINT"
                print "&END MOTION"
            }
            if (calc_type == "VIBRATIONAL_ANALYSIS") {
                print "&VIBRATIONAL_ANALYSIS"
                print "  DX 0.01"
                print "  NPROC_REP 4" 
                print "  TC_PRESSURE 101325" 
                print "  TC_TEMPERATURE 298.15"
                print "  THERMOCHEMISTRY"
                print "  INTENSITIES F"
                print "  FULLY_PERIODIC T"
                print "  &PRINT"
                print "    &MOLDEN_VIB" 
                print "    &END MOLDEN_VIB"
                print "  &END PRINT"
                print "&END VIBRATIONAL_ANALYSIS"
            }
        }
        ' "$CIF2CELL_IN"

    mv $CIF2CELL_IN ${CIF_FILE%.*}.$CALC_TYPE.inp

}

# convert QE-in file to cif file
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

# prepare QE-in files for conv test with the provided QE-in file
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

# analysis QE-out file for scf conv
out_scf () {

    OUT_FILE="$1"
    
    echo "Cycle,Energy,Accuracy" > "${OUT_FILE%.*}_scf.csv"
    paste -d ' ' \
        <(grep  "iteration #" "$OUT_FILE" | awk -F'[#]' '{printf "%d\n", $2}') \
        <(grep  "total energy              =" "$OUT_FILE" | awk -F'[=]' '{printf "%.10f\n", $2}') \
        <(grep  "estimated scf accuracy" "$OUT_FILE" | awk -F'[<]' '{printf "%.6e\n", $2}') \
    | awk '{printf "%d,%.10f,%.6e\n", $1, $2, $3}' >> "${OUT_FILE%.*}_scf.csv"

}

# analysis QE-out file for relax/vc-relax conv
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

# generate the cif file through cell lattice and atomic coordinates from publications 
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

    CIF_GEN_FILE="${1:-structure}"

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

# ==============================================================================
# Main Execution Logic
# ==============================================================================

# Check for help flag immediately
if [[ $# -eq 0 || $1 == "-h" || $1 == "--help" ]]; then
    print_help
    exit 0
fi

# Validate file input or prompt if invalid
READ_FILE="$1"

# Check if file is missing
if [[ ! -f "$READ_FILE" ]]; then
    
    # logic: Only prompt to generate if the missing file ends in .cif
    if [[ "$READ_FILE" == *.cif ]]; then
        echo "CIF file '$READ_FILE' was not found."
        read -p "Do you want to generate it manually? [Y/n]: " cif_gen_choice
        cif_gen_choice=${cif_gen_choice:-Y}
        
        if [[ "$cif_gen_choice" =~ ^[Yy]$ ]]; then
            cif_gen "$READ_FILE"
            exit 0
        fi
    fi

    # If it's not a .cif file, or the user said "No", exit with error
    echo "Error: Input file '$READ_FILE' not found."
    exit 1
fi

# Check if file exists but is not readable (permissions error)
if [[ ! -r "$READ_FILE" ]]; then
    echo "Error: Input file '$READ_FILE' exists but is not readable."
    exit 1
fi

# Detect file type
EXT="${READ_FILE##*.}"
case "$EXT" in
    cif) 
        echo "Successfully loaded CIF file: $READ_FILE"
        CIF_FILE=$READ_FILE
        APP_FORMAT="${2:-qe}"
        CALC_TYPE="${3:-scf}"

        case "$APP_FORMAT" in
            qe) 
                echo "Starting to prepare the $CALC_TYPE input files for Quantum ESPRESSO program"
                cif2qe $CIF_FILE $CALC_TYPE
                ;;
            vasp) 
                echo "Starting to prepare the $CALC_TYPE input files for VASP program"
                cif2vasp $CIF_FILE $CALC_TYPE
                ;;
            cp2k) 
                echo "Starting to prepare the $CALC_TYPE input files for CP2K program"
                cif2cp2k $CIF_FILE $CALC_TYPE
                ;;
            *) 
                echo "Unknown application format: $APP_FORMAT"
                exit 1
                ;;
        esac
        ;;
    in)  
        echo "Successfully loaded Quantum ESPRESSO IN file: $READ_FILE"; 
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
        echo "Successfully loaded Quantum ESPRESSO OUT file: $READ_FILE"; 
        OUT_FILE=$READ_FILE 
        OUT_ANALYSIS="${2:-relax}"

        case "$OUT_ANALYSIS" in
            scf) 
                echo "Analysis the Quantum ESPRESSO scf output files"
                out_scf $OUT_FILE
                ;;
            relax|vc-relax) 
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
