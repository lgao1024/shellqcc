#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage:" >&2
  echo "  Single file: $0 input_file [output_file]" >&2
  echo "    - If input suffix is .cif/.CIF => CIF -> TXT" >&2
  echo "    - Otherwise => TXT -> CIF (fallback behavior)" >&2
  echo "  Multi file:  $0 input1 input2 [input3 ...]" >&2
  echo "    - Type is determined by input1 suffix" >&2
  echo "    - If input1 is .cif: convert all *.cif inputs to *.txt" >&2
  echo "    - If input1 is not .cif: generate all TXT combinations to CIFs" >&2
  exit 2
}

[[ $# -lt 1 ]] && usage

# ---------------- helpers (shared) ----------------

clean_line() {
  local s="$1"
  s="${s%%#*}"
  s="$(printf '%s' "$s" | awk '{$1=$1; print}')"
  printf '%s' "$s"
}

num_to_float() {
  local s="$1"
  awk -v x="$s" '
    function trim(t){ gsub(/^[ \t\r\n]+|[ \t\r\n]+$/, "", t); return t }
    BEGIN{
      x=trim(x)
      sub(/\([0-9]+\)$/, "", x)   # 0.123(4) -> 0.123

      if (x ~ /^[+-]?[0-9]+\/[0-9]+$/) {
        n=x; sub(/\/.*/, "", n)
        d=x; sub(/.*\//, "", d)
        if (d==0) exit 3
        printf "%.12g", (n+0.0)/(d+0.0)
        exit 0
      }
      printf "%.12g", (x+0.0)     # awk parses .5 too
    }
  '
}

is_number_like() {
  local s="$1"
  awk -v x="$s" '
    function trim(t){ gsub(/^[ \t\r\n]+|[ \t\r\n]+$/, "", t); return t }
    BEGIN{
      x=trim(x)
      sub(/\([0-9]+\)$/, "", x)
      if (x ~ /^[+-]?[0-9]+\/[0-9]+$/) { print "Y"; exit }
      if (x ~ /^[+-]?([0-9]+(\.[0-9]*)?|\.[0-9]+)$/) { print "Y"; exit }
      print "N"
    }
  ' | grep -q '^Y$'
}

fmt6() {
  local s="$1"
  awk -v x="$s" 'BEGIN{ printf "%.6f", (x+0.0) }'
}

# Parse first line: sg a b c alpha beta gamma
# - sg may be quoted with '...' or "..." (to allow spaces, e.g. 'P 21/c')
# Output via echo: sg|a|b|c|alpha|beta|gamma
parse_text_header() {
  local line="$1"
  awk -v L="$line" '
    function trim(t){ gsub(/^[ \t\r\n]+|[ \t\r\n]+$/, "", t); return t }
    function strip_uncert(s){ sub(/\([0-9]+\)$/, "", s); return s }
    function tofloat(s,   n,d){
      s=trim(s); s=strip_uncert(s)
      if (s ~ /^[+-]?[0-9]+\/[0-9]+$/) {
        n=s; sub(/\/.*/, "", n)
        d=s; sub(/.*\//, "", d)
        return (n+0.0)/(d+0.0)
      }
      return (s+0.0)
    }
    BEGIN{
      L=trim(L)
      # extract sg
      if (L ~ /^'\''/) {
        sub(/^'\''/, "", L)
        sg=L
        sub(/'\''.*/, "", sg)
        rest=L
        sub(/^[^'\'']*'\''[ \t]*/, "", rest)
      } else if (L ~ /^"/) {
        sub(/^"/, "", L)
        sg=L
        sub(/".*/, "", sg)
        rest=L
        sub(/^[^"]*"[ \t]*/, "", rest)
      } else {
        sg=L
        sub(/[ \t].*$/, "", sg)
        rest=L
        sub(/^[^ \t]+[ \t]*/, "", rest)
      }

      n=split(rest, t, /[ \t]+/)
      if (n < 6) {
        print "ERR"
        exit
      }
      a=tofloat(t[1]); b=tofloat(t[2]); c=tofloat(t[3])
      al=tofloat(t[4]); be=tofloat(t[5]); ga=tofloat(t[6])

      print sg "|" a "|" b "|" c "|" al "|" be "|" ga
    }
  '
}

# ---------------- TXT -> CIF ----------------

txt2cif() {
  local in="$1"
  local out="$2"

  local sg="" a="" b="" c="" alpha="" beta="" gamma=""

  declare -A elem_count
  declare -a A_label A_elem A_x A_y A_z A_occ
  local atom_n=0
  local first_done=0

  while IFS= read -r raw || [[ -n "$raw" ]]; do
    local line
    line="$(clean_line "$raw")"
    [[ -z "$line" ]] && continue

    if [[ $first_done -eq 0 ]]; then
      local parsed
      parsed="$(parse_text_header "$line")"
      if [[ "$parsed" == "ERR" || -z "$parsed" ]]; then
        echo "ERROR: first line must be: sg a b c alpha beta gamma" >&2
        echo "Got: $line" >&2
        exit 2
      fi
      IFS="|" read -r sg a b c alpha beta gamma <<< "$parsed"
      first_done=1
      continue
    fi

    IFS=$' \t' read -r -a tok <<< "$line"
    if [[ ${#tok[@]} -lt 4 ]]; then
      echo "ERROR: atom line has too few columns: $line" >&2
      exit 2
    fi

    local label="" elem="" x="" y="" z="" occ=""

    # If 2nd token is numeric => element x y z [occ]
    if is_number_like "${tok[1]}"; then
      elem="${tok[0]}"
      x="$(num_to_float "${tok[1]}")"
      y="$(num_to_float "${tok[2]}")"
      z="$(num_to_float "${tok[3]}")"
      if [[ ${#tok[@]} -ge 5 ]]; then
        occ="$(num_to_float "${tok[4]}")"
      else
        occ="1"
      fi

      elem_count["$elem"]=$(( ${elem_count["$elem"]:-0} + 1 ))
      label="${elem}${elem_count["$elem"]}"
    else
      # label element x y z [occ]
      if [[ ${#tok[@]} -lt 5 ]]; then
        echo "ERROR: expected at least 5 fields: label element x y z [occ]" >&2
        echo "Got: $line" >&2
        exit 2
      fi
      label="${tok[0]}"
      elem="${tok[1]}"
      x="$(num_to_float "${tok[2]}")"
      y="$(num_to_float "${tok[3]}")"
      z="$(num_to_float "${tok[4]}")"
      if [[ ${#tok[@]} -ge 6 ]]; then
        occ="$(num_to_float "${tok[5]}")"
      else
        occ="1"
      fi
    fi

    A_label[$atom_n]="$label"
    A_elem[$atom_n]="$elem"
    A_x[$atom_n]="$x"
    A_y[$atom_n]="$y"
    A_z[$atom_n]="$z"
    A_occ[$atom_n]="$occ"
    atom_n=$((atom_n + 1))
  done < "$in"

  if [[ $first_done -eq 0 ]]; then
    echo "ERROR: no valid first line found in $in" >&2
    exit 2
  fi

  local data_name="${out##*/}"
  data_name="${data_name%.*}"
  data_name="${data_name// /_}"

  {
    echo "data_${data_name}"
    echo
    echo "_symmetry_space_group_name_H-M  '$(printf "%s" "$sg")'"
    echo
    echo "_cell_length_a    $(fmt6 "$a")"
    echo "_cell_length_b    $(fmt6 "$b")"
    echo "_cell_length_c    $(fmt6 "$c")"
    echo "_cell_angle_alpha $(fmt6 "$alpha")"
    echo "_cell_angle_beta  $(fmt6 "$beta")"
    echo "_cell_angle_gamma $(fmt6 "$gamma")"
    echo
    echo "loop_"
    echo "_atom_site_label"
    echo "_atom_site_type_symbol"
    echo "_atom_site_fract_x"
    echo "_atom_site_fract_y"
    echo "_atom_site_fract_z"
    echo "_atom_site_occupancy"

    for ((i=0; i<atom_n; i++)); do
      printf "%s %s %s %s %s %s\n" \
        "${A_label[$i]}" \
        "${A_elem[$i]}" \
        "$(fmt6 "${A_x[$i]}")" \
        "$(fmt6 "${A_y[$i]}")" \
        "$(fmt6 "${A_z[$i]}")" \
        "$(fmt6 "${A_occ[$i]}")"
    done
    echo
  } > "$out"

  echo "Wrote: $out"
}

# ---------------- CIF -> TXT ----------------
# Tries to read:
#  - space group: _symmetry_space_group_name_H-M or _space_group_name_H-M_alt
#  - cell: _cell_length_a/b/c, _cell_angle_alpha/beta/gamma
#  - atom loop: a loop_ containing _atom_site_fract_x/_y/_z and type_symbol (optional) + label/occupancy (optional)
cif2txt() {
  local in="$1"
  local out="$2"

  awk '
    function trim(s){ gsub(/^[ \t\r\n]+|[ \t\r\n]+$/, "", s); return s }
    function strip_uncert(s){ sub(/\([0-9]+\)$/, "", s); return s }
    function normnum(s,   n,d){
      s=trim(s); s=strip_uncert(s)
      if (s ~ /^[+-]?[0-9]+\/[0-9]+$/) {
        n=s; sub(/\/.*/, "", n)
        d=s; sub(/.*\//, "", d)
        if (d == 0) return 0
        return (n+0.0)/(d+0.0)
      }
      return (s+0.0)
    }
    function unquote(s){
      s=trim(s)
      if (s ~ /^'\''.*'\''$/) { sub(/^'\''/, "", s); sub(/'\''$/, "", s); return s }
      if (s ~ /^".*"$/) { sub(/^"/, "", s); sub(/"$/, "", s); return s }
      return s
    }
    function fmt6(x){ return sprintf("%.6f", x+0.0) }

    BEGIN{
      sg=""
      a=b=c=alpha=beta=gamma=""
      in_loop=0; reading_data=0
      ncol=0
      idx_label=idx_type=idx_x=idx_y=idx_z=idx_occ=0
      have_atoms=0
      use_label_col=0
      # store atom lines for printing after header
      n_atoms=0
    }

    # remove comments starting with #
    {
      line=$0
      sub(/#.*/, "", line)
      line=trim(line)
      if (line == "") next
    }

    # capture space group
    tolower(line) ~ /^_symmetry_space_group_name_h-m[ \t]/ {
      val=line
      sub(/^[^ \t]+[ \t]+/, "", val)
      sg=unquote(val)
      next
    }
    tolower(line) ~ /^_space_group_name_h-m_alt[ \t]/ && sg=="" {
      val=line
      sub(/^[^ \t]+[ \t]+/, "", val)
      sg=unquote(val)
      next
    }

    # capture cell params
    tolower(line) ~ /^_cell_length_a[ \t]/ { val=line; sub(/^[^ \t]+[ \t]+/, "", val); a=normnum(val); next }
    tolower(line) ~ /^_cell_length_b[ \t]/ { val=line; sub(/^[^ \t]+[ \t]+/, "", val); b=normnum(val); next }
    tolower(line) ~ /^_cell_length_c[ \t]/ { val=line; sub(/^[^ \t]+[ \t]+/, "", val); c=normnum(val); next }
    tolower(line) ~ /^_cell_angle_alpha[ \t]/ { val=line; sub(/^[^ \t]+[ \t]+/, "", val); alpha=normnum(val); next }
    tolower(line) ~ /^_cell_angle_beta[ \t]/ { val=line; sub(/^[^ \t]+[ \t]+/, "", val); beta=normnum(val); next }
    tolower(line) ~ /^_cell_angle_gamma[ \t]/ { val=line; sub(/^[^ \t]+[ \t]+/, "", val); gamma=normnum(val); next }

    # loop parsing
    tolower(line) == "loop_" {
      in_loop=1
      reading_data=0
      ncol=0
      idx_label=idx_type=idx_x=idx_y=idx_z=idx_occ=0
      next
    }

    # collect loop headers
    in_loop && substr(line,1,1)=="_" && reading_data==0 {
      ncol++
      col[ncol]=tolower(line)
      if (col[ncol] ~ /^_atom_site_label$/) idx_label=ncol
      if (col[ncol] ~ /^_atom_site_type_symbol$/) idx_type=ncol
      if (col[ncol] ~ /^_atom_site_fract_x$/) idx_x=ncol
      if (col[ncol] ~ /^_atom_site_fract_y$/) idx_y=ncol
      if (col[ncol] ~ /^_atom_site_fract_z$/) idx_z=ncol
      if (col[ncol] ~ /^_atom_site_occupancy$/) idx_occ=ncol
      next
    }

    # first data line after headers (if this loop is atom loop)
    in_loop && reading_data==0 && substr(line,1,1)!="_" {
      # determine if this loop is the atom loop we want
      if (idx_x>0 && idx_y>0 && idx_z>0) {
        reading_data=1
        use_label_col = (idx_label>0)
        # fall through to parse this line as data
      } else {
        # not an atom loop; stop tracking
        in_loop=0
        reading_data=0
        next
      }
    }

    # read atom data lines
    in_loop && reading_data==1 {
      # stop conditions: new loop, new tag, new data block
      if (tolower(line)=="loop_" || substr(line,1,1)=="_" || tolower(substr(line,1,5))=="data_") {
        in_loop=0
        reading_data=0
        next
      }
      # split on whitespace
      nf=split(line, f, /[ \t]+/)
      if (nf < ncol) next  # skip malformed line

      # extract fields
      label = (idx_label>0 ? f[idx_label] : "")
      elem  = (idx_type>0 ? f[idx_type] : "")
      x = normnum(f[idx_x]); y = normnum(f[idx_y]); z = normnum(f[idx_z])

      occ = "NA"
      if (idx_occ>0) {
        occ_raw=f[idx_occ]
        if (occ_raw!="." && occ_raw!="?") occ = normnum(occ_raw)
      }

      # build output line:
      # if label exists => label elem x y z [occ?]
      # else => elem x y z [occ?]
      if (elem=="" && label!="") {
        # derive element from label (strip trailing digits)
        elem=label
        sub(/[0-9]+$/, "", elem)
      }
      if (elem=="") elem="X"

      lineout=""
      if (use_label_col) {
        lineout = label " " elem " " fmt6(x) " " fmt6(y) " " fmt6(z)
      } else {
        lineout = elem " " fmt6(x) " " fmt6(y) " " fmt6(z)
      }
      if (occ!="NA") lineout = lineout " " fmt6(occ)

      atoms[++n_atoms]=lineout
      have_atoms=1
      next
    }

    END{
      if (sg=="") sg="P1"
      if (a==""||b==""||c==""||alpha==""||beta==""||gamma=="") {
        # minimal fallback if missing
        if (a=="") a=1; if (b=="") b=1; if (c=="") c=1
        if (alpha=="") alpha=90; if (beta=="") beta=90; if (gamma=="") gamma=90
      }

      # quote sg if it contains spaces
      if (sg ~ /[ \t]/) sg_out="'"'"'" sg "'"'"'"
      else sg_out=sg

      print sg_out, fmt6(a), fmt6(b), fmt6(c), fmt6(alpha), fmt6(beta), fmt6(gamma)

      if (have_atoms) {
        for (i=1; i<=n_atoms; i++) print atoms[i]
      }
    }
  ' "$in" > "$out"

  echo "Wrote: $out"
}

# ---------------- batch helpers ----------------

get_type_prefix() {
  local file="$1"
  local base="${file##*/}"
  base="${base%.*}"
  local prefix
  prefix="$(printf '%s' "$base" | sed -E 's/[0-9].*$//')"
  if [[ -z "$prefix" ]]; then
    prefix="$base"
  fi
  printf '%s' "$prefix"
}

get_first_nonempty_clean_line() {
  local file="$1"
  while IFS= read -r raw || [[ -n "$raw" ]]; do
    local line
    line="$(clean_line "$raw")"
    [[ -z "$line" ]] && continue
    printf '%s\n' "$line"
    return 0
  done < "$file"
  return 1
}

collect_atom_lines_from_txt() {
  local file="$1"
  local first_done=0
  while IFS= read -r raw || [[ -n "$raw" ]]; do
    local line
    line="$(clean_line "$raw")"
    [[ -z "$line" ]] && continue
    if [[ $first_done -eq 0 ]]; then
      first_done=1
      continue
    fi
    printf '%s\n' "$line"
  done < "$file"
}

format_and_symmetrize_cifs() {
  local src_dir="$1"
  local formatted_dir="${src_dir}/formatted"
  local symm_dir="${src_dir}/symm"
  local formatted_count=0
  local symm_count=0

  local have_ase=1
  local have_c2x=1
  command -v ase >/dev/null 2>&1 || have_ase=0
  command -v c2x >/dev/null 2>&1 || have_c2x=0

  if [[ $have_ase -eq 0 && $have_c2x -eq 0 ]]; then
    echo "WARNING: neither 'ase' nor 'c2x' found, skip formatting/symmetrization."
    return 0
  fi
  if [[ $have_ase -eq 0 ]]; then
    echo "WARNING: 'ase' not found, skip formatted CIF generation."
  fi
  if [[ $have_c2x -eq 0 ]]; then
    echo "WARNING: 'c2x' not found, skip symmetrized CIF generation."
  fi

  mkdir -p "$formatted_dir" "$symm_dir"

  local -a cif_files=()
  local f
  shopt -s nullglob
  cif_files=("$src_dir"/*.cif)
  shopt -u nullglob

  if [[ ${#cif_files[@]} -eq 0 ]]; then
    echo "No CIF files found in $src_dir for formatting."
    return 0
  fi

  for f in "${cif_files[@]}"; do
    local base out_fmt out_symm symm_in
    base="${f##*/}"
    out_fmt="$formatted_dir/$base"
    out_symm="$symm_dir/${base%.*}_symm.cif"

    if [[ $have_ase -eq 1 ]]; then
      if ase convert -i cif "$f" -o cif "$out_fmt"; then
        formatted_count=$((formatted_count + 1))
      else
        echo "WARNING: ase convert failed for $f"
      fi
    fi

    if [[ $have_c2x -eq 1 ]]; then
      if [[ $have_ase -eq 1 && -f "$out_fmt" ]]; then
        symm_in="$out_fmt"
      else
        symm_in="$f"
      fi
      if c2x -r --cif --sym "$symm_in" "$out_symm"; then
        symm_count=$((symm_count + 1))
      else
        echo "WARNING: c2x symmetrization failed for $symm_in"
      fi
    fi
  done

  echo "Formatted CIFs generated: $formatted_count -> $formatted_dir"
  echo "Symmetrized CIFs generated: $symm_count -> $symm_dir"
}

batch_cif_to_txt() {
  local -a all_inputs=("$@")
  local -a cif_inputs=()
  local f
  for f in "${all_inputs[@]}"; do
    [[ -f "$f" ]] || { echo "ERROR: input not found: $f" >&2; exit 2; }
    local lf="${f##*/}"
    lf="${lf,,}"
    if [[ "$lf" == *.cif ]]; then
      cif_inputs+=("$f")
    fi
  done

  [[ ${#cif_inputs[@]} -gt 0 ]] || { echo "ERROR: first input is CIF mode, but no .cif files provided." >&2; exit 2; }

  for f in "${cif_inputs[@]}"; do
    cif2txt "$f" "${f%.*}.txt"
  done
}

batch_txt_to_cif_combinations() {
  local -a all_inputs=("$@")
  local -a txt_inputs=()
  local f
  for f in "${all_inputs[@]}"; do
    [[ -f "$f" ]] || { echo "ERROR: input not found: $f" >&2; exit 2; }
    local lf="${f##*/}"
    lf="${lf,,}"
    if [[ "$lf" == *.txt ]]; then
      txt_inputs+=("$f")
    fi
  done

  [[ ${#txt_inputs[@]} -gt 0 ]] || { echo "ERROR: first input is TXT mode, but no .txt files provided." >&2; exit 2; }

  local first_header
  first_header="$(get_first_nonempty_clean_line "${txt_inputs[0]}")" || {
    echo "ERROR: no valid header line found in ${txt_inputs[0]}" >&2
    exit 2
  }
  local parsed_header
  parsed_header="$(parse_text_header "$first_header")"
  if [[ "$parsed_header" == "ERR" || -z "$parsed_header" ]]; then
    echo "ERROR: first TXT file header is invalid: ${txt_inputs[0]}" >&2
    echo "Header: $first_header" >&2
    exit 2
  fi
  local sg a b c alpha beta gamma
  IFS="|" read -r sg a b c alpha beta gamma <<< "$parsed_header"

  declare -A PREFIX_FILES
  declare -A PREFIX_COUNTS
  declare -A PREFIX_ORDER_SET
  local -a PREFIX_ORDER=()
  local p

  for f in "${txt_inputs[@]}"; do
    p="$(get_type_prefix "$f")"
    if [[ -z "${PREFIX_ORDER_SET[$p]:-}" ]]; then
      PREFIX_ORDER_SET["$p"]=1
      PREFIX_ORDER+=("$p")
      PREFIX_FILES["$p"]="$f"
    else
      PREFIX_FILES["$p"]+=$'\x1f'"$f"
    fi
  done

  local prefix
  for prefix in "${PREFIX_ORDER[@]}"; do
    IFS=$'\x1f' read -r -a _files <<< "${PREFIX_FILES[$prefix]}"
    PREFIX_COUNTS["$prefix"]="${#_files[@]}"
  done

  # Pre-parse each TXT once:
  # output: elem x y z occ
  # labels in source TXT are intentionally ignored for speed and uniqueness.
  declare -A FILE_ATOM_TMP
  local -a TMP_FILES=()
  for f in "${txt_inputs[@]}"; do
    local atom_tmp
    atom_tmp="$(mktemp)"
    TMP_FILES+=("$atom_tmp")
    FILE_ATOM_TMP["$f"]="$atom_tmp"

    local first_done=0
    while IFS= read -r raw || [[ -n "$raw" ]]; do
      local line
      line="$(clean_line "$raw")"
      [[ -z "$line" ]] && continue

      if [[ $first_done -eq 0 ]]; then
        first_done=1
        continue
      fi

      IFS=$' \t' read -r -a tok <<< "$line"
      [[ ${#tok[@]} -lt 4 ]] && continue

      local elem x y z occ
      if is_number_like "${tok[1]}"; then
        # element x y z [occ]
        elem="${tok[0]}"
        x="$(num_to_float "${tok[1]}")"
        y="$(num_to_float "${tok[2]}")"
        z="$(num_to_float "${tok[3]}")"
        if [[ ${#tok[@]} -ge 5 ]]; then
          occ="$(num_to_float "${tok[4]}")"
        else
          occ="1"
        fi
      else
        # label element x y z [occ] -> ignore source label
        [[ ${#tok[@]} -lt 5 ]] && continue
        elem="${tok[1]}"
        x="$(num_to_float "${tok[2]}")"
        y="$(num_to_float "${tok[3]}")"
        z="$(num_to_float "${tok[4]}")"
        if [[ ${#tok[@]} -ge 6 ]]; then
          occ="$(num_to_float "${tok[5]}")"
        else
          occ="1"
        fi
      fi
      printf "%s %s %s %s %s\n" "$elem" "$x" "$y" "$z" "$occ" >> "$atom_tmp"
    done < "$f"
  done

  local out_dir="txt_combo_cifs"
  mkdir -p "$out_dir"

  local total_combinations=1
  for prefix in "${PREFIX_ORDER[@]}"; do
    total_combinations=$(( total_combinations * PREFIX_COUNTS["$prefix"] ))
  done

  local comb_idx
  for ((comb_idx=0; comb_idx<total_combinations; comb_idx++)); do
    local remaining_idx="$comb_idx"
    local -a chosen_files=()
    local -a chosen_bases=()

    for prefix in "${PREFIX_ORDER[@]}"; do
      IFS=$'\x1f' read -r -a _files <<< "${PREFIX_FILES[$prefix]}"
      local count="${#_files[@]}"
      local file_idx=$(( remaining_idx % count ))
      remaining_idx=$(( remaining_idx / count ))
      local selected="${_files[$file_idx]}"
      chosen_files+=("$selected")
      local selected_base="${selected##*/}"
      chosen_bases+=("${selected_base%.*}")
    done

    local joined_name
    local IFS='-'
    joined_name="${chosen_bases[*]}"
    local out_path="$out_dir/$joined_name.cif"

    {
      echo "data_${joined_name}"
      echo
      echo "_symmetry_space_group_name_H-M  '$(printf "%s" "$sg")'"
      echo
      echo "_cell_length_a    $(fmt6 "$a")"
      echo "_cell_length_b    $(fmt6 "$b")"
      echo "_cell_length_c    $(fmt6 "$c")"
      echo "_cell_angle_alpha $(fmt6 "$alpha")"
      echo "_cell_angle_beta  $(fmt6 "$beta")"
      echo "_cell_angle_gamma $(fmt6 "$gamma")"
      echo
      echo "loop_"
      echo "_atom_site_label"
      echo "_atom_site_type_symbol"
      echo "_atom_site_fract_x"
      echo "_atom_site_fract_y"
      echo "_atom_site_fract_z"
      echo "_atom_site_occupancy"
    } > "$out_path"

    local -a atom_chunks=()
    local cf
    for cf in "${chosen_files[@]}"; do
      atom_chunks+=("${FILE_ATOM_TMP[$cf]}")
    done

    awk '
      {
        elem=$1; x=$2; y=$3; z=$4; occ=$5
        n[elem]++
        label=elem n[elem]
        printf "%s %s %.6f %.6f %.6f %.6f\n", label, elem, x+0.0, y+0.0, z+0.0, occ+0.0
      }
    ' "${atom_chunks[@]}" >> "$out_path"
  done

  rm -f "${TMP_FILES[@]}"
  echo "Generated $total_combinations CIF files in: $out_dir"
  format_and_symmetrize_cifs "$out_dir"
}

# ---------------- dispatch by suffix & argc ----------------

if [[ $# -eq 1 ]]; then
  in="$1"
  [[ -f "$in" ]] || { echo "ERROR: input not found: $in" >&2; exit 2; }
  lower_in="${in##*/}"
  lower_in="${lower_in,,}"
  if [[ "$lower_in" == *.cif ]]; then
    cif2txt "$in" "${in%.*}.txt"
  else
    txt2cif "$in" "${in%.*}.cif"
  fi
elif [[ $# -eq 2 ]]; then
  first="${1##*/}"
  first="${first,,}"
  second="${2##*/}"
  second="${second,,}"
  # If two args share the same source suffix kind, treat as two-input batch mode.
  if { [[ "$first" == *.cif ]] && [[ "$second" == *.cif ]]; } || \
     { [[ "$first" != *.cif ]] && [[ "$second" == *.txt ]]; }; then
    if [[ "$first" == *.cif ]]; then
      batch_cif_to_txt "$@"
    else
      batch_txt_to_cif_combinations "$@"
    fi
  else
    in="$1"
    out="$2"
    [[ -f "$in" ]] || { echo "ERROR: input not found: $in" >&2; exit 2; }
    lower_in="${in##*/}"
    lower_in="${lower_in,,}"
    if [[ "$lower_in" == *.cif ]]; then
      cif2txt "$in" "$out"
    else
      txt2cif "$in" "$out"
    fi
  fi
else
  first="${1##*/}"
  first="${first,,}"
  if [[ "$first" == *.cif ]]; then
    batch_cif_to_txt "$@"
  else
    batch_txt_to_cif_combinations "$@"
  fi
fi
