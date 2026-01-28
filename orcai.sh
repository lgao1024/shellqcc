#!/bin/bash
# orcai.sh — Automate ORCA calculations
# Usage: ./orcai.sh <template_file>

# ==============================================================================
# Global Configuration
# ==============================================================================
ORCA_BIN='/opt/orca-6.1.1/orca'
TEMPLATE_FILE="$1"

# ==============================================================================
# Function Definitions
# ==============================================================================

function xtb_prepare() {
    local template="$1"
    echo "[Mode] Running XTB Preparation with template: $template"

    # 1. Convert .mol to .xyz (if missing)
    for mol_file in *.mol; do
        [[ -f "$mol_file" ]] || continue
        
        base_name="${mol_file%.*}"
        xyz_file="${base_name}.xyz"

        if [[ ! -f "$xyz_file" ]]; then
            echo "[Convert] Generating $xyz_file from $mol_file..."
            PYTHONWARNINGS="ignore" ase convert "$mol_file" "$xyz_file"
        fi
    done

    # 2. Loop through XYZ files to prepare and run
    for xyz_file in *.xyz; do
        [[ -f "$xyz_file" ]] || continue
        
        base_name="${xyz_file%.*}.xtb"
        mol_file="${xyz_file%.*}.mol"  # Define the potential mol file
        inp_file="${base_name}.inp"
        out_file="${base_name}.out"
        xtbxyz_file="${base_name}.xyz"
        work_dir="${base_name}"

        echo "[Prepare] Processing $xyz_file..."

        # Generate Input: specific sed replacement for XTB template
        sed "/\* XYZFILE/s|$| $xyz_file|" "$template" > "$inp_file"

        # Create folder
        if [[ ! -d "$work_dir" ]]; then
            mkdir -p "$work_dir"
        fi
        
        # Move files
        mv "$xyz_file" "$inp_file" "$work_dir/"
        if [[ -f "$mol_file" ]]; then
            mv "$mol_file" "$work_dir/"
        fi

        # Run ORCA
        (
            cd "$work_dir" || exit
            echo "[Running] Starting ORCA (XTB) for $base_name..."
            "$ORCA_BIN" "$inp_file" |tee "$out_file"
            if [[ -f "$xtbxyz_file" ]]; then
                cp "$xtbxyz_file" "../"
            fi
        )
    done
}

function input_prepare() {
    local template="$1"
    echo "[Mode] Running Standard Input Preparation with template: $template"
    
    for xyz_file in *.xyz; do
        [[ -f "$xyz_file" ]] || continue
        base_name="${xyz_file%%.*}.${template%.*}"
        inp_file="${base_name}.inp"
        out_file="${base_name}.out"
        sed "/\* XYZFILE/s|$| ${xyz_file%.*}.xyz|" "$template" > "$base_name.inp"
        # "$ORCA_BIN" "$inp_file" |tee "$out_file"
    done
}

# ==============================================================================
# Main Execution Logic
# ==============================================================================

# 1. Check if argument exists
if [[ -z "$TEMPLATE_FILE" ]]; then
    echo "Error: No template file provided."
    echo "Usage: $0 <template_file>"
    exit 1
fi

# 2. Check if file exists (Handle missing xtb_tmp specifically)
if [[ ! -f "$TEMPLATE_FILE" ]]; then
    
    # Specific logic for missing 'xtb_tmp'
    if [[ "$TEMPLATE_FILE" == "xtb.tmp" ]]; then
        echo "Template 'xtb.tmp' not found. Generating default xtb.tmp file."

        cat > xtb.tmp <<EOF
%maxcore 1024
%pal nprocs 4 end
! Opt XTB SlowConv DefGrid3 TightSCF MiniPrint

%geom
   MaxIter 1024
end

* XYZFILE 0 1
EOF
        echo "Created 'xtb.tmp'."

        # Ask user how to proceed
        echo "Select action:"
        echo "  1) Continue execution (Run calculations now)"
        echo "  2) Exit (Just generate file)"
        read -p "Choice [1]: " action_choice
        action_choice=${action_choice:-1}

        if [[ "$action_choice" != "1" ]]; then
            echo "Exiting."
            exit 0
        fi
    else
        echo "Error: Template file '$TEMPLATE_FILE' not found."
        exit 1
    fi
fi

# 3. Route Logic based on Filename
if [[ "$TEMPLATE_FILE" == "xtb.tmp" ]]; then
    xtb_prepare "$TEMPLATE_FILE"
else
    input_prepare "$TEMPLATE_FILE"
fi

echo "All tasks finished."