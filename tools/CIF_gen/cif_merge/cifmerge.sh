#!/bin/bash

# Configuration
TEMPLATE="cif.template"
TEMP_ATOMS="temp_combined_atoms.txt"
MERGED_DIR="./merged_cifs"
FORMATTED_DIR="./merged_cifs_formatted"
merged_count=0   # Counter for successfully created merged CIFs
formatted_count=0 # Counter for successfully formatted CIFs


# --------------------------
# Step 0: Create output folders
# --------------------------
mkdir -p "$MERGED_DIR" || { echo "Error: Failed to create $MERGED_DIR!"; exit 1; }
mkdir -p "$FORMATTED_DIR" || { echo "Error: Failed to create $FORMATTED_DIR!"; exit 1; }
mkdir -p "$FORMATTED_DIR/symm" || { echo "Error: Failed to create $FORMATTED_DIR/symm!"; exit 1; }
echo "Merged files will be saved to: $MERGED_DIR"
echo "Formatted files will be saved to: $FORMATTED_DIR"


# --------------------------
# Step 1: Detect structure types
# --------------------------
get_prefix() {
    local filename=$(basename "$1" .cif)
    echo "$filename" | sed -E 's/[0-9].*//'
}

declare -A PREFIXES
for cif in *.cif; do
    [ -f "$cif" ] || continue
    prefix=$(get_prefix "$cif")
    [ -n "$prefix" ] && PREFIXES["$prefix"]=1
done

SORTED_PREFIXES=($(printf "%s\n" "${!PREFIXES[@]}" | sort))

if [ ${#SORTED_PREFIXES[@]} -lt 2 ]; then
    echo "Error: Need at least 2 structure types (e.g., A*.cif and B*.cif)."
    exit 1
fi

echo "Detected structure types (sorted): ${SORTED_PREFIXES[*]}"


# --------------------------
# Step 2: Collect files for each type
# --------------------------
declare -A TYPE_FILES
declare -A TYPE_COUNTS

for prefix in "${SORTED_PREFIXES[@]}"; do
    files=()
    for f in "${prefix}"*.cif; do
        [ -f "$f" ] && files+=("$f")
    done
    
    if [ ${#files[@]} -eq 0 ]; then
        echo "Error: No files found for type '$prefix'."
        exit 1
    fi
    
    TYPE_FILES["$prefix"]="${files[*]}"
    TYPE_COUNTS["$prefix"]=${#files[@]}
    echo "Type $prefix: ${#files[@]} files (${files[*]})"
done


# --------------------------
# Step 3: Generate, count, and format all combinations
# --------------------------
total_combinations=1
for prefix in "${SORTED_PREFIXES[@]}"; do
    total_combinations=$((total_combinations * TYPE_COUNTS["$prefix"]))
done

echo "Total combinations to process: $total_combinations"

echo -e "\n=========================================================================================="

# Generate each combination
for ((comb_idx=0; comb_idx < total_combinations; comb_idx++)); do
    declare -A CURRENT_FILES=()
    declare -A CURRENT_IDS=()
    
    remaining_idx=$comb_idx
    
    for prefix in "${SORTED_PREFIXES[@]}"; do
        files_str="${TYPE_FILES[$prefix]}"
        files=($files_str)
        count=${#files[@]}
        
        file_idx=$((remaining_idx % count))
        remaining_idx=$((remaining_idx / count))
        
        current_file="${files[$file_idx]}"
        CURRENT_FILES["$prefix"]="$current_file"
        
        base_name=$(basename "$current_file" .cif)
        id="${base_name#$prefix}"
        CURRENT_IDS["$prefix"]="$id"
    done


    # Build output filename with underscores (e.g., A1_B2.cif)
    output_filename=""
    for prefix in "${SORTED_PREFIXES[@]}"; do
        output_filename+="${prefix}${CURRENT_IDS[$prefix]}_"
    done
    output_filename="${output_filename%_}.cif"


    # Combine atomic data and create merged file
    > "$TEMP_ATOMS"
    for prefix in "${SORTED_PREFIXES[@]}"; do
        file="${CURRENT_FILES[$prefix]}"
        echo "Extracting atoms from $file..."
        awk '/_atom_site_occupancy/ {flag=1; next} /loop_/ {flag=0} flag' "$file" | \
        awk '{print $2, $3, $4, $5, $8}' >> "$TEMP_ATOMS"
    done


    # Save merged file and verify creation
    merged_path="$MERGED_DIR/$output_filename"
    sed -e "/_atom_site_occupancy/ r $TEMP_ATOMS" "$TEMPLATE" > "$merged_path"
    
    if [ -f "$merged_path" ]; then  # Only count if file exists
        ((merged_count++))
        echo "==> Generated merged file: $merged_path (Count: $merged_count)"
    else
        echo "Warning: Failed to create merged file $merged_path"
        continue  # Skip formatting if merged file failed
    fi


    # Reformat with ASE and verify success
    formatted_path="$FORMATTED_DIR/$output_filename"
    formatted_symm_path="$FORMATTED_DIR/symm/${output_filename%.*}_symm.cif"
    if ase convert -i cif "$merged_path" -o cif "$formatted_path"; then
        ((formatted_count++))
        echo "==> Generated formatted file: $formatted_path (Count: $formatted_count)"
        c2x -r --cif --sym "$formatted_path" "$formatted_symm_path"
        echo "==> Generated formatted file: $formatted_symm_path (Count: $formatted_count)"
        echo -e "==========================================================================================\n"
    else
        echo "Warning: Failed to format $merged_path (not counted)"
    fi
done


# Cleanup and final report
rm "$TEMP_ATOMS"

echo -e "\n============================================="
echo "                FINAL RESULTS                "
echo "============================================="
echo "  Merged CIFs created:       $merged_count"
echo "  Stored in:                 $MERGED_DIR"
echo "---------------------------------------------"
echo "  Formatted CIFs generated:  $formatted_count"
echo "  Stored in:                 $FORMATTED_DIR"
echo "  Stored in (symmetrized):   $FORMATTED_DIR/symm"
echo "============================================="
echo -e "Done! Check the directories above for output files.\n"


