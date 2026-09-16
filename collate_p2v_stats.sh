#!/usr/bin/env bash

set -euo pipefail

root="${1:-output}"

if [[ ! -d "$root" ]]; then
    echo "No such directory: $root" >&2
    exit 1
fi

mapfile -d '' stats_files < <(
    find "$root" -type f \
        -name 'p2v_stats_*modes*trials*over*seeing*arcsec*.dat' \
        -print0 | sort -z -V
)

if (( ${#stats_files[@]} == 0 )); then
    echo "No p2v stats files found under $root" >&2
    exit 1
fi

generated_list="$(mktemp)"
trap 'rm -f "$generated_list"' EXIT

awk -v generated_list="$generated_list" '
function reset_file()
{
    source_dir = ""
    aperture = ""
    modes = ""
    trials = ""
    oversize = ""
    seeing = ""
    p2v_mean = ""
    p2v_stddev = ""
    p2v_max = ""
    p2pdiff_mean = ""
    p2pdiff_stddev = ""
    p2pdiff_max = ""
    skip_file = 0
}

function parse_filename( basename )
{
    basename = current_file
    sub(/^.*\//, "", basename)

    source_dir = current_file
    sub(/\/[^\/]*$/, "", source_dir)

    # The aperture may contain underscores, and the filename may have a
    # run-specific suffix after the seeing value.
    if (basename !~ /^p2v_stats_.+_[0-9]+modes_[0-9]+trials_[0-9]+over_seeing[0-9.]+arcsec(_.*)?\.dat$/) {
        print "Skipping unrecognized stats filename: " current_file > "/dev/stderr"
        skip_file = 1
        return
    }

    aperture = basename
    sub(/^p2v_stats_/, "", aperture)
    sub(/_[0-9]+modes_.*/, "", aperture)

    modes = basename
    sub(/^p2v_stats_/, "", modes)
    sub("^" aperture "_", "", modes)
    sub(/modes_.*/, "", modes)

    trials = basename
    sub("^p2v_stats_" aperture "_[0-9]+modes_", "", trials)
    sub(/trials_.*/, "", trials)

    oversize = basename
    sub("^p2v_stats_" aperture "_[0-9]+modes_[0-9]+trials_", "", oversize)
    sub(/over_seeing.*/, "", oversize)

    seeing = basename
    sub(/^.*_seeing/, "", seeing)
    sub(/arcsec(_.*)?\.dat$/, "", seeing)
}

function write_header_once( out_file )
{
    if (!(out_file in wrote_header)) {
        print "# modes p2v_5sigma p2v_max p2pdiff_5sigma p2pdiff_max" > out_file
        wrote_header[out_file] = 1
        print out_file >> generated_list
    }
}

function flush_file( out_file, p2v_5sigma, p2pdiff_5sigma )
{
    if (current_file == "" || skip_file) {
        return
    }

    if (p2v_mean == "" || p2v_stddev == "" || p2v_max == "" ||
        p2pdiff_mean == "" || p2pdiff_stddev == "" || p2pdiff_max == "") {
        print "Skipping incomplete stats file: " current_file > "/dev/stderr"
        return
    }

    out_file = source_dir "/p2v_stats_" aperture "_seeing" seeing "arcsec.dat"
    p2v_5sigma = p2v_mean + 5 * p2v_stddev
    p2pdiff_5sigma = p2pdiff_mean + 5 * p2pdiff_stddev

    write_header_once(out_file)

    printf "%d %.10g %.10g %.10g %.10g\n", \
        modes + 0, p2v_5sigma, p2v_max, p2pdiff_5sigma, p2pdiff_max >> out_file
}

FNR == 1 {
    flush_file()
    current_file = FILENAME
    reset_file()
    parse_filename()
}

$1 == "p2v" {
    p2v_mean = $2
    p2v_stddev = (NF >= 5 ? $4 : $3)
    # Current files have: metric mean rms stddev max.  Accept the older
    # four-field form, where the third value is the dispersion, as well.
    p2v_max = (NF >= 5 ? $5 : $4)
}

$1 == "p2pdiff" {
    p2pdiff_mean = $2
    p2pdiff_stddev = (NF >= 5 ? $4 : $3)
    p2pdiff_max = (NF >= 5 ? $5 : $4)
}

END {
    flush_file()
}
' "${stats_files[@]}"

if [[ ! -s "$generated_list" ]]; then
    echo "No CSV files generated." >&2
    exit 1
fi

echo "Wrote:"
sort -V "$generated_list"
