#!/bin/bash

# ==========================================
# Help Information and Input Validation
# ==========================================
if [ "$#" -eq 0 ]; then
    echo "Usage:"
    echo "  1. Split file: $0 <merged_cif_file.cif>"
    echo "  2. Merge files: $0 <file1.cif> <file2.cif> ... <fileN.cif>"
    exit 1
fi

# ==========================================
# Function 1: Split a single merged CIF file
# ==========================================
if [ "$#" -eq 1 ]; then
    INPUT_FILE="$1"
    
    if [ ! -f "$INPUT_FILE" ]; then
        echo "Error: File '$INPUT_FILE' not found"
        exit 1
    fi

    # Get the base filename without extension and create output directory
    BASENAME="${INPUT_FILE%.*}"
    OUT_DIR="${BASENAME}_splited"
    mkdir -p "$OUT_DIR"

    echo "Starting to split CIF files into directory: $OUT_DIR/"

    # Use awk to read and output block by block
    awk -v out_dir="$OUT_DIR" '
    BEGIN {
        counter = 0; has_data = 0; main_buffer = ""; temp_buffer = ""; ccdc_num = ""
    }
    /^data_/ {
        if (has_data) {
            filename = out_dir "/" ((ccdc_num != "") ? ccdc_num : "structure" counter) ".cif"
            print main_buffer > filename
            close(filename)
            print "Generated: " filename
        }
        counter++; ccdc_num = ""; has_data = 1
        
        if (temp_buffer != "") main_buffer = temp_buffer "\n" $0
        else main_buffer = $0
        temp_buffer = ""
        next
    }
    {
        if (!has_data) {
            if (temp_buffer == "") temp_buffer = $0
            else temp_buffer = temp_buffer "\n" $0
            next
        }
        if ($0 ~ /^[[:space:]]*#/ || $0 ~ /^[[:space:]]*$/) {
            if (temp_buffer == "") temp_buffer = $0
            else temp_buffer = temp_buffer "\n" $0
        } else {
            if (temp_buffer != "") {
                main_buffer = main_buffer "\n" temp_buffer
                temp_buffer = ""
            }
            main_buffer = main_buffer "\n" $0
            
            if ($1 == "_database_code_depnum_ccdc_archive") {
                match($0, /[0-9]+/)
                if (RSTART > 0) ccdc_num = substr($0, RSTART, RLENGTH)
            }
        }
    }
    END {
        if (has_data) {
            if (temp_buffer != "") main_buffer = main_buffer "\n" temp_buffer
            filename = out_dir "/" ((ccdc_num != "") ? ccdc_num : "structure" counter) ".cif"
            print main_buffer > filename
            close(filename)
            print "Generated: " filename
        }
    }
    ' "$INPUT_FILE"
    
    echo "Split complete!"
    exit 0
fi

# ==========================================
# Function 2: Merge multiple CIF files
# ==========================================
files=("$@")
# Check if all input files exist
for f in "${files[@]}"; do
    if [ ! -f "$f" ]; then
        echo "Error: File '$f' not found"
        exit 1
    fi
done

# Initialize default order array
current_order=()
for (( i=1; i<=${#files[@]}; i++ )); do
    current_order+=("$i")
done

# Function to parse user input (comma-separated and hyphen-ranges)
parse_input_to_order() {
    local input="${1//,/ }" # Replace commas with spaces
    local new_order=()
    for item in $input; do
        # Match hyphen range, e.g., 5-6 or 3-1
        if [[ "$item" =~ ^[0-9]+-[0-9]+$ ]]; then
            start=${item%-*}
            end=${item#*-}
            if [ "$start" -le "$end" ]; then
                for (( j=start; j<=end; j++ )); do
                    new_order+=("$j")
                done
            else
                for (( j=start; j>=end; j-- )); do
                    new_order+=("$j")
                done
            fi
        # Match pure numbers
        elif [[ "$item" =~ ^[0-9]+$ ]]; then
            new_order+=("$item")
        fi
    done
    echo "${new_order[@]}"
}

# Interactive sorting loop
while true; do
    echo -e "\nCurrent files to be merged and their order:"
    for i in "${!current_order[@]}"; do
        idx=${current_order[$i]}
        real_idx=$((idx - 1))
        echo "  $((i + 1)): ${files[$real_idx]}"
    done

    echo -e "\nYou can enter the desired merge order (e.g., 1, 3, 5-6, 2)."
    read -p "Press Enter to keep the current order and merge directly: " user_input

    # If empty input, break the loop and start merging
    if [[ -z "$user_input" ]]; then
        break
    fi

    # Parse user input into a new order array
    parsed_output=$(parse_input_to_order "$user_input")
    read -ra temp_order <<< "$parsed_output"

    # Validate the input indices
    valid=true
    for idx in "${temp_order[@]}"; do
        if [ "$idx" -lt 1 ] || [ "$idx" -gt "${#files[@]}" ]; then
            echo "Warning: Index $idx is out of range! Please try again."
            valid=false
            break
        fi
    done

    if [ "$valid" = true ]; then
        # Preview the new order
        echo -e "\n---- Preview of updated file order ----"
        for i in "${!temp_order[@]}"; do
            idx=${temp_order[$i]}
            real_idx=$((idx - 1))
            echo "  $((i + 1)): ${files[$real_idx]}"
        done
        
        read -p "Confirm to merge with this order? (y/n, Enter to confirm): " confirm
        if [[ -z "$confirm" || "$confirm" =~ ^[Yy]$ ]]; then
            current_order=("${temp_order[@]}")
            break
        fi
    fi
done

# Execute merge
MERGED_OUTPUT="merged_output_$(date +%Y%m%d_%H%M%S).cif"
echo -e "\nStarting merge..."
> "$MERGED_OUTPUT" # Clear or create the output file

for i in "${!current_order[@]}"; do
    idx=${current_order[$i]}
    real_idx=$((idx - 1))
    cat "${files[$real_idx]}" >> "$MERGED_OUTPUT"
    # Ensure at least one newline separator between files
    echo -e "\n" >> "$MERGED_OUTPUT"
done

echo "Merge complete! Successfully output to file: $MERGED_OUTPUT"