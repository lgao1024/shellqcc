#!/usr/bin/env bash
set -euo pipefail

# Auto-generate a styled gnuplot script from tabular data and export PDF/PNG.
#
# Usage:
#   ./qeplot.sh data.csv
#   ./qeplot.sh data.csv xy
#   ./qeplot.sh data.csv xyy
#   ./qeplot.sh data.csv xyyy
#
# Mode semantics:
#   xy      : one panel, all Y columns share the same Y axis
#   xyy+    : stacked multiplot panels, one independent Y axis per Y column
#
# Notes:
# - Relative paths only
# - First row: legends
# - Second row: units if mostly non-numeric
# - Generated outputs:
#     <stem>_auto_plot.gnuplot
#     <stem>.pdf
#     <stem>.png

if [ "$#" -lt 1 ] || [ "$#" -gt 2 ]; then
  echo "Usage: $0 <datafile> [xy|xyy|xyyy|...]" >&2
  exit 1
fi

DATAFILE=$1
MODE=${2:-xy}

if [ ! -f "$DATAFILE" ]; then
  echo "Error: file not found: $DATAFILE" >&2
  exit 1
fi

if ! command -v gnuplot >/dev/null 2>&1; then
  echo "Error: gnuplot is not installed or not in PATH." >&2
  exit 1
fi

case "$MODE" in
  xy) PANEL_MODE=0 ;;
  xy*)
    case "$MODE" in
      x) echo "Error: mode must include at least one y" >&2; exit 1 ;;
      *[!xy]*) echo "Error: mode must look like xy, xyy, xyyy, ..." >&2; exit 1 ;;
      *) PANEL_MODE=1 ;;
    esac
    ;;
  *)
    echo "Error: mode must look like xy, xyy, xyyy, ..." >&2
    exit 1
    ;;
esac

BASENAME=$(basename -- "$DATAFILE")
STEM=${BASENAME%.*}
GNUPLOT_SCRIPT="${STEM}_auto_plot.gnuplot"
PDF_OUT="${STEM}.pdf"
PNG_OUT="${STEM}.png"

detect_delim() {
  local file sample comma_count tab_count semicolon_count pipe_count
  file=$1
  sample=$(head -n 2 "$file" || true)

  comma_count=$(printf '%s\n' "$sample" | awk -F',' '{s+=NF-1} END{print s+0}')
  tab_count=$(printf '%s\n' "$sample" | awk -F'\t' '{s+=NF-1} END{print s+0}')
  semicolon_count=$(printf '%s\n' "$sample" | awk -F';' '{s+=NF-1} END{print s+0}')
  pipe_count=$(printf '%s\n' "$sample" | awk -F'|' '{s+=NF-1} END{print s+0}')

  if [ "$tab_count" -gt "$comma_count" ] && [ "$tab_count" -ge "$semicolon_count" ] && [ "$tab_count" -ge "$pipe_count" ] && [ "$tab_count" -gt 0 ]; then
    printf '\t'
  elif [ "$comma_count" -ge "$semicolon_count" ] && [ "$comma_count" -ge "$pipe_count" ] && [ "$comma_count" -gt 0 ]; then
    printf ','
  elif [ "$semicolon_count" -ge "$pipe_count" ] && [ "$semicolon_count" -gt 0 ]; then
    printf ';'
  elif [ "$pipe_count" -gt 0 ]; then
    printf '|'
  else
    printf 'whitespace'
  fi
}

DELIM=$(detect_delim "$DATAFILE")

count_columns() {
  local file delim
  file=$1
  delim=$2
  if [ "$delim" = "whitespace" ]; then
    awk 'NF>0 {print NF; exit}' "$file"
  else
    awk -v FS="$delim" 'NF>0 {print NF; exit}' "$file"
  fi
}

NCOLS=$(count_columns "$DATAFILE" "$DELIM")
if [ -z "$NCOLS" ] || [ "$NCOLS" -lt 2 ]; then
  echo "Error: could not detect at least two columns in $DATAFILE" >&2
  exit 1
fi

NY=$((NCOLS - 1))

SECOND_IS_UNITS=$(awk -v mode="$DELIM" '
function isnum(s) {
  gsub(/^ +| +$/, "", s)
  return (s ~ /^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/)
}
NR==2 {
  n=0; non=0
  if (mode=="whitespace") {
    for (i=1; i<=NF; i++) {
      n++
      if (!isnum($i)) non++
    }
  } else {
    split($0, a, mode)
    for (i=1; i<=length(a); i++) {
      n++
      gsub(/^ +| +$/, "", a[i])
      if (!isnum(a[i])) non++
    }
  }
  if (n>0 && non >= int((n+1)/2)) print 1; else print 0
  exit
}' "$DATAFILE")

if [ -z "$SECOND_IS_UNITS" ]; then
  SECOND_IS_UNITS=0
fi

DATA_START_LINE=2
if [ "$SECOND_IS_UNITS" -eq 1 ]; then
  DATA_START_LINE=3
fi

get_field() {
  local line_no col_no
  line_no=$1
  col_no=$2
  if [ "$DELIM" = "whitespace" ]; then
    awk -v ln="$line_no" -v c="$col_no" 'NR==ln {print $c; exit}' "$DATAFILE"
  else
    awk -v FS="$DELIM" -v ln="$line_no" -v c="$col_no" 'NR==ln {print $c; exit}' "$DATAFILE"
  fi
}

trim_quotes_spaces() {
  sed -e 's/^ *//' -e 's/ *$//' -e 's/^"//' -e 's/"$//'
}

make_label() {
  local idx header unit
  idx=$1
  header=$(get_field 1 "$idx" | trim_quotes_spaces)
  if [ -z "$header" ]; then
    if [ "$idx" -eq 1 ]; then
      header="x"
    else
      header="y"
    fi
  fi

  if [ "$SECOND_IS_UNITS" -eq 1 ]; then
    unit=$(get_field 2 "$idx" | trim_quotes_spaces)
    if [ -n "$unit" ] && [ "$unit" != "-" ] && [ "$unit" != "NA" ] && [ "$unit" != "na" ]; then
      printf '%s (%s)' "$header" "$unit"
      return
    fi
  fi
  printf '%s' "$header"
}

escape_gp() {
  printf '%s' "$1" | sed "s/'/''/g"
}

X_LABEL=$(make_label 1)

build_common_header() {
  cat <<EOF
reset
set encoding utf8
set datafile commentschars "#"
EOF

  if [ "$DELIM" != "whitespace" ]; then
    printf "set datafile separator '%s'\n" "$DELIM"
  fi

  cat <<'EOF'

set border lw 1.2
set tics nomirror out scale 0.75
set xtics font ",10"
set ytics font ",10"
set mxtics 2
set mytics 2
set grid xtics ytics mxtics mytics lw 0.6 lc rgb "#d9d9d9"
set key outside right top vertical box opaque spacing 1.15 width 0.5
set key font ",10"

set style line 1  lc rgb "#1f77b4" lw 2.0 pt 7  ps 0.9
set style line 2  lc rgb "#d62728" lw 2.0 pt 5  ps 0.9
set style line 3  lc rgb "#2ca02c" lw 2.0 pt 9  ps 0.9
set style line 4  lc rgb "#9467bd" lw 2.0 pt 13 ps 0.9
set style line 5  lc rgb "#ff7f0e" lw 2.0 pt 11 ps 0.9
set style line 6  lc rgb "#17becf" lw 2.0 pt 3  ps 0.9
set style line 7  lc rgb "#8c564b" lw 2.0 pt 4  ps 0.9
set style line 8  lc rgb "#e377c2" lw 2.0 pt 6  ps 0.9
set style line 9  lc rgb "#7f7f7f" lw 2.0 pt 8  ps 0.9
set style line 10 lc rgb "#bcbd22" lw 2.0 pt 10 ps 0.9

# Optional manual controls:
# set xrange [xmin:xmax]
# set yrange [ymin:ymax]
# set logscale y
# set size ratio 0.6

EOF
}

write_single_panel_body() {
  local plot_cmd c
  plot_cmd=''
  c=2
  while [ "$c" -le "$NCOLS" ]; do
    if [ -z "$plot_cmd" ]; then
      plot_cmd="'${DATAFILE}' using 1:${c} every ::$((DATA_START_LINE-1)) with linespoints ls ${c-1} title columnhead(${c})"
    else
      plot_cmd="${plot_cmd},\\
    '${DATAFILE}' using 1:${c} every ::$((DATA_START_LINE-1)) with linespoints ls ${c-1} title columnhead(${c})"
    fi
    c=$((c+1))
  done

  cat <<EOF
set xlabel '$(escape_gp "$X_LABEL")' font ",12"
set ylabel 'y' font ",12"
set title '$(escape_gp "$STEM")' font ",13"

plot ${plot_cmd}

EOF
}

write_stacked_panels_body() {
  local c idx label bottom screen_height top bottom_pos
  idx=1
  screen_height=$(awk -v n="$NY" 'BEGIN { printf "%.6f", 1.0/n }')

  cat <<EOF
unset key
set multiplot layout ${NY},1 rowsfirst title '$(escape_gp "$STEM")' font ",13"
EOF

  c=2
  while [ "$c" -le "$NCOLS" ]; do
    label=$(make_label "$c")
    if [ "$c" -lt "$NCOLS" ]; then
      cat <<EOF
unset xlabel
set format x ""
EOF
    else
      cat <<EOF
set xlabel '$(escape_gp "$X_LABEL")' font ",12"
set format x
EOF
    fi

    cat <<EOF
set ylabel '$(escape_gp "$label")' font ",11"
plot '${DATAFILE}' using 1:${c} every ::$((DATA_START_LINE-1)) with linespoints ls ${idx} title columnhead(${c})

EOF
    idx=$((idx+1))
    c=$((c+1))
  done

  cat <<'EOF'
unset multiplot

EOF
}

write_gnuplot_file() {
  local terminal_block
  terminal_block=$1
  {
    build_common_header
    printf "%s\n" "$terminal_block"
    if [ "$PANEL_MODE" -eq 0 ]; then
      write_single_panel_body
    else
      write_stacked_panels_body
    fi
  } > "$GNUPLOT_SCRIPT"
}

PDF_HEIGHT=$(awk -v n="$NY" 'BEGIN { h = 2.0 + 2.2*n; if (h < 4.6) h = 4.6; printf "%.2f", h }')
PNG_HEIGHT=$(awk -v n="$NY" 'BEGIN { h = 350 + 320*n; if (h < 1200) h = 1200; printf "%d", h }')

write_gnuplot_file "set terminal pdfcairo enhanced color font \",11\" size 6.8in,${PDF_HEIGHT}in
set output '${PDF_OUT}'"

gnuplot "$GNUPLOT_SCRIPT"

{
  build_common_header
  printf "set terminal pngcairo enhanced truecolor font \",11\" size 1800,%s\n" "$PNG_HEIGHT"
  printf "set output '%s'\n\n" "$PNG_OUT"
  if [ "$PANEL_MODE" -eq 0 ]; then
    write_single_panel_body
  else
    write_stacked_panels_body
  fi
} > "$GNUPLOT_SCRIPT"

gnuplot "$GNUPLOT_SCRIPT"

echo "Generated gnuplot script: $GNUPLOT_SCRIPT"
echo "Generated PDF:            $PDF_OUT"
echo "Generated PNG:            $PNG_OUT"
