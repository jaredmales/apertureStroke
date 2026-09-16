#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"

n_trials="${1:-10000}"
n_threads="${2:-0}"
output_root="${3:-output/gmagaox_stroke_exploration_hybrid}"

segment_modes="${SEGMENT_MODES:-/home/jrmales/Source/mxWork/GMT/DMs/stroke/svd/p2v_25m_segment_zernike_pinv/gmt_segment_modes_182x182.fits}"

seeings=(0.5 0.64 0.79)
outer_scales=(25 50 100 1000 0)

# These primary-mode cutoffs are the nine low Zernikes plus the rectangular
# modified-Fourier sequence. For each linear DOF N, mxlib uses
# 4*floor(N/2)*(floor(N/2)+1) Fourier planes.
hybrid_linear_dofs=(1 2 4 8 16 32)
hybrid_zernikes=9
hybrid_cutoff_values=(0)
for linear_dof in "${hybrid_linear_dofs[@]}"; do
    half_dof=$((linear_dof / 2))
    fourier_modes=$((4 * half_dof * (half_dof + 1)))
    hybrid_cutoff_values+=( $((hybrid_zernikes + fourier_modes)) )
done
hybrid_cutoffs="$(IFS=,; printf '%s' "${hybrid_cutoff_values[*]}")"
hybrid_final_cutoff="${hybrid_cutoff_values[${#hybrid_cutoff_values[@]} - 1]}"

zernike_6p5_cutoffs="0,5,10,20,40,60,80,100"

if [[ ! -x ./apertureStroke || ! -x ./zernikeTheory ]]; then
    echo "Build apertureStroke and zernikeTheory first: make" >&2
    exit 1
fi

if [[ ! -f "$segment_modes" ]]; then
    echo "Segment PTT cube not found: $segment_modes" >&2
    echo "Set SEGMENT_MODES=/path/to/gmt_segment_modes_182x182.fits" >&2
    exit 1
fi

scalar_tag()
{
    local value="$1"
    value="${value//./p}"
    value="${value//-/m}"
    printf '%s' "$value"
}

outer_scale_tag()
{
    local outer_scale="$1"
    if [[ "$outer_scale" == "0" ]]; then
        printf '%s' "L0inf"
    else
        printf 'L0%sm' "$(scalar_tag "$outer_scale")"
    fi
}

run_theory()
{
    local config="$1"
    local case_name="$2"
    local n_modes="$3"
    local seeing="$4"
    local seeing_tag
    seeing_tag="$(scalar_tag "$seeing")"
    local theory_dir="$output_root/theory/$case_name/seeing$seeing_tag"
    local theory_label="${case_name}_seeing${seeing_tag}"
    local theory_file="$theory_dir/zernike_theory_${theory_label}.dat"

    if [[ "${FORCE:-0}" != "1" && -s "$theory_file" ]]; then
        echo "Skipping existing theory $theory_file"
        return
    fi

    ./zernikeTheory -c "$config" \
        --basis.modes="$n_modes" \
        --atmosphere.seeing="$seeing" \
        --output.directory="$theory_dir" \
        --output.label="$theory_label"
}

run_case()
{
    local config="$1"
    local case_name="$2"
    local label="$3"
    local cutoffs="$4"
    local final_cutoff="$5"
    local prefix_file="${6:-}"

    for outer_scale in "${outer_scales[@]}"; do
        local outer_tag
        outer_tag="$(outer_scale_tag "$outer_scale")"

        for seeing in "${seeings[@]}"; do
            local seeing_tag
            seeing_tag="$(scalar_tag "$seeing")"
            local output_dir="$output_root/$case_name/$outer_tag/seeing$seeing_tag"

            shopt -s nullglob
            local existing=("$output_dir"/p2v_stats_"$label"_"$final_cutoff"modes_*.dat)
            shopt -u nullglob
            if [[ "${FORCE:-0}" != "1" && ${#existing[@]} -gt 0 ]]; then
                echo "Skipping existing $output_dir"
                continue
            fi

            cmd=(./apertureStroke -c "$config"
                 "--simulation.trials=$n_trials"
                 "--simulation.threads=$n_threads"
                 "--atmosphere.seeing=$seeing"
                 "--atmosphere.outerScale=$outer_scale"
                 "--basis.cutoffs=$cutoffs"
                 "--output.directory=$output_dir")
            if [[ -n "$prefix_file" ]]; then
                cmd+=("--basis.prefixFile=$prefix_file" "--basis.prefixName=segmentPTT")
            fi

            printf 'Running:'
            printf ' %q' "${cmd[@]}"
            printf '\n'
            "${cmd[@]}"
        done
    done
}

echo "Hybrid primary cutoffs (linear DOF ${hybrid_linear_dofs[*]}): $hybrid_cutoffs"

# The 6.5 m references remain Zernike bases. Their theory is generated here
# because it remains directly applicable to those runs.
for seeing in "${seeings[@]}"; do
    run_theory exploration_circle_6p5m.conf circle_6p5m 100 "$seeing"
done

run_case exploration_circle_6p5m.conf circle_6p5m circle_6p5m \
    "$zernike_6p5_cutoffs" 100
run_case exploration_magaox_6p5m.conf magaox_6p5m magaox_6p5m \
    "$zernike_6p5_cutoffs" 100
run_case exploration_circle_25m_hybrid.conf circle_25m_hybrid_z2to10_fourierRectN32 \
    circle_25m_hybrid_z2to10_fourierRectN32 "$hybrid_cutoffs" "$hybrid_final_cutoff"
run_case exploration_gmt_25m_hybrid.conf gmt_25m_hybrid_z2to10_fourierRectN32 \
    gmt_25m_hybrid_z2to10_fourierRectN32 "$hybrid_cutoffs" "$hybrid_final_cutoff"
run_case exploration_gmt_25m_segmentPTT_hybridN32.conf gmt_25m_segmentPTT_hybrid_z2to10_fourierRectN32 \
    gmt_25m_segmentPTT_hybrid_z2to10_fourierRectN32 "$hybrid_cutoffs" "$hybrid_final_cutoff" "$segment_modes"
