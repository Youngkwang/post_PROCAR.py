#!/usr/bin/env python3

"""VASP PROCAR post-processing script.

Parses PROCAR and OUTCAR files to generate band structure outputs
in Igor Pro .itx, matplotlib, and/or CSV formats.

Modernized from post_PROCAR.py (2015-2022).
"""

import math
import os
import sys

import vasp_io as vasp
import output_writers as writers
from config import get_config


# ---------------------------------------------------------------------------
# Parsing functions
# ---------------------------------------------------------------------------

def parse_kvector(outcar_lines, num_kpoints, starting_kpoint, outcar_file):
    """Parse k-vector cumulative distances from OUTCAR.

    Returns a list of length num_kpoints with cumulative k-space distances.
    The first (starting_kpoint - 1) entries remain 0.

    Physics-critical: preserves exact k-vector distance calculation.
    """
    two_pi = 2.0 * math.pi

    kvec_line = outcar_lines.index(
        vasp.get_keyword_from("2pi/SCALE", outcar_file)
    )
    kvec_loc_bef = kvec_line + 1

    parts_bef = outcar_lines[kvec_loc_bef].split()
    kvec_X_bef = float(parts_bef[0]) * two_pi
    kvec_Y_bef = float(parts_bef[1]) * two_pi
    kvec_Z_bef = float(parts_bef[2]) * two_pi
    kvec_K_bef = 0.0

    kvec_mat = [0.0] * num_kpoints

    for i in range(starting_kpoint, num_kpoints + 1):
        j = kvec_line + i
        parts = outcar_lines[j].split()

        kvec_X = float(parts[0]) * two_pi
        kvec_Y = float(parts[1]) * two_pi
        kvec_Z = float(parts[2]) * two_pi

        kvec_K = kvec_K_bef + math.sqrt(
            (kvec_X - kvec_X_bef) ** 2
            + (kvec_Y - kvec_Y_bef) ** 2
            + (kvec_Z - kvec_Z_bef) ** 2
        )
        kvec_mat[i - 1] = kvec_K

        kvec_X_bef = kvec_X
        kvec_Y_bef = kvec_Y
        kvec_Z_bef = kvec_Z
        kvec_K_bef = kvec_K

    return kvec_mat


def parse_procar(procar_lines, num_kpoints, num_bands, num_ions,
                 starting_band, ending_band, starting_ion, ending_ion,
                 starting_kpoint, check_soc, ref_energy):
    """Parse eigenvalues, occupations, and orbital projections from PROCAR.

    Returns (eigen_val_mat, elec_occ_mat, pro_band_mat).

    Physics-critical: PROCAR line-index formulas for SOC vs non-SOC are
    preserved exactly from the original script.
    """
    eigen_val_mat = [[0.0] * num_kpoints for _ in range(num_bands)]
    elec_occ_mat = [[0.0] * num_kpoints for _ in range(num_bands)]
    # pro_band_mat[orbital_type][ion][band][kpoint]
    pro_band_mat = [
        [[[0.0] * num_kpoints for _ in range(num_bands)]
         for _ in range(num_ions + 1)]
        for _ in range(4)
    ]

    for m in range(starting_band, ending_band + 1):
        print("Band number %03d is being generated" % m)

        for o in range(starting_ion, ending_ion + 2):

            for k in range(starting_kpoint, num_kpoints + 1):
                # Physics-critical index formulas — preserved exactly
                if check_soc == "F":
                    h = (4 + (k - 1) * (num_bands * (num_ions + 5) + 3)
                         + (2 + (m - 1) * (num_ions + 5)))
                elif check_soc == "T":
                    h = (4 + (k - 1) * (num_bands * 4 * (num_ions + 2) + 3)
                         + (2 + (m - 1) * 4 * (num_ions + 2)))
                else:
                    raise ValueError(
                        "Unexpected LSORBIT value: %s" % check_soc)

                eigen_parts = procar_lines[h - 1].split()
                eigen_VAL = float(eigen_parts[4])
                eigen_val_mat[m - 1][k - 1] = eigen_VAL - ref_energy
                elec_occ_mat[m - 1][k - 1] = float(eigen_parts[7])

                if check_soc == "F":
                    l = (4 + (k - 1) * (num_bands * (num_ions + 5) + 3)
                         + (4 + (m - 1) * (num_ions + 5)) + o)
                elif check_soc == "T":
                    l = (4 + (k - 1) * (num_bands * 4 * (num_ions + 2) + 3)
                         + (4 + (m - 1) * 4 * (num_ions + 2)) + o)

                proj_parts = procar_lines[l - 1].split()
                pro_band_mat[0][o - 1][m - 1][k - 1] = float(proj_parts[1])
                pro_band_mat[1][o - 1][m - 1][k - 1] = float(proj_parts[2])
                pro_band_mat[2][o - 1][m - 1][k - 1] = float(proj_parts[3])
                pro_band_mat[3][o - 1][m - 1][k - 1] = float(proj_parts[4])

    return eigen_val_mat, elec_occ_mat, pro_band_mat


def aggregate_projections_by_element(pro_band_mat, element_type, element_num,
                                     orbital_types, num_kpoints, num_bands,
                                     starting_band, ending_band,
                                     starting_kpoint, ending_kpoint):
    """Aggregate per-ion projections into per-element projections.

    Returns pro_band_elements_mat[orbital][element][band][kpoint].

    Physics-critical: element aggregation loop preserved exactly.
    """
    pro_band_elements_mat = [
        [[[0.0] * num_kpoints for _ in range(num_bands)]
         for _ in range(len(element_type))]
        for _ in range(4)
    ]

    for t in range(len(orbital_types)):
        element_num_start = 0
        for z in range(len(element_type)):
            element_num_end = element_num_start + int(element_num[z])
            for e in range(element_num_start, element_num_end):
                for y in range(starting_band - 1, ending_band):
                    for x in range(starting_kpoint - 1, ending_kpoint):
                        pro_band_elements_mat[t][z][y][x] += (
                            pro_band_mat[t][e][y][x]
                        )
            element_num_start = element_num_end

    return pro_band_elements_mat


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    cfg = get_config()

    # Validate input files exist
    for key in ("procar_file", "outcar_file"):
        if not os.path.isfile(cfg[key]):
            print("Error: file not found: %s" % cfg[key], file=sys.stderr)
            sys.exit(1)

    # Read input files with proper resource management
    with open(cfg["procar_file"], "r") as f:
        procar_lines = f.readlines()
    with open(cfg["outcar_file"], "r") as f:
        outcar_lines = f.readlines()

    # Parse header info from PROCAR
    header = procar_lines[1].split()
    num_ions = int(header[11])
    num_kpoints = int(header[3])
    num_bands = int(header[7])

    # Reference energy: auto-detect from OUTCAR E-fermi, or use provided value
    if cfg["ref_energy"] == "auto":
        fermi_line = vasp.get_keyword_from("E-fermi", cfg["outcar_file"])
        if fermi_line is None:
            print("Error: could not find E-fermi in OUTCAR", file=sys.stderr)
            sys.exit(1)
        ref_energy = float(fermi_line.split()[2])
        print("Auto-detected E-fermi: %.5f eV" % ref_energy)
    else:
        ref_energy = float(cfg["ref_energy"])

    # Fill in defaults for ranges that depend on parsed values
    starting_kpoint = cfg["starting_kpoint"]
    ending_kpoint = (cfg["ending_kpoint"]
                     if cfg["ending_kpoint"] is not None
                     else num_kpoints)
    starting_band = cfg["starting_band"]
    ending_band = (cfg["ending_band"]
                   if cfg["ending_band"] is not None
                   else num_bands)
    starting_ion = cfg["starting_ion"]
    ending_ion = (cfg["ending_ion"]
                  if cfg["ending_ion"] is not None
                  else num_ions)

    # Update cfg with resolved values for writers
    cfg["ending_kpoint"] = ending_kpoint
    cfg["ending_band"] = ending_band
    cfg["ending_ion"] = ending_ion

    # Check SOC setting
    soc_line = vasp.get_keyword_from("LSORBIT", cfg["outcar_file"])
    if soc_line is None:
        print("Error: could not find LSORBIT in OUTCAR", file=sys.stderr)
        sys.exit(1)
    check_soc = soc_line.split()[2]

    # Parse k-vector distances
    print("K-vector is being generated")
    kvec_mat = parse_kvector(
        outcar_lines, num_kpoints, starting_kpoint, cfg["outcar_file"])

    # Parse PROCAR data
    eigen_val_mat, elec_occ_mat, pro_band_mat = parse_procar(
        procar_lines, num_kpoints, num_bands, num_ions,
        starting_band, ending_band, starting_ion, ending_ion,
        starting_kpoint, check_soc, ref_energy,
    )

    # Element info for per-element projections
    contcar = cfg["contcar_file"]
    if not os.path.isfile(contcar):
        print("Error: file not found: %s" % contcar, file=sys.stderr)
        sys.exit(1)
    element_type = vasp.get_elements_type(contcar)
    element_num = vasp.get_elements_num(contcar)
    element_info = vasp.get_elements_info(contcar)

    print("num_IONS: %d" % num_ions)
    print("num_KPOINTS: %d" % num_kpoints)
    print("num_BANDS: %d" % num_bands)
    print("ref_E: %.4f eV" % ref_energy)
    print(element_info)

    # Aggregate projections by element
    orbital_types = cfg["orbital_types"]
    pro_band_elements_mat = aggregate_projections_by_element(
        pro_band_mat, element_type, element_num, orbital_types,
        num_kpoints, num_bands, starting_band, ending_band,
        starting_kpoint, ending_kpoint,
    )

    # Generate outputs
    output_formats = cfg["output_formats"]

    if "igor" in output_formats:
        print("Writing Igor Pro .itx files...")
        writers.write_igor_band_structure(kvec_mat, eigen_val_mat, cfg)
        writers.write_igor_occ_band(
            kvec_mat, eigen_val_mat, elec_occ_mat, cfg)
        writers.write_igor_pband_per_ion(
            kvec_mat, eigen_val_mat, pro_band_mat, cfg)
        writers.write_igor_pband_per_element(
            kvec_mat, eigen_val_mat, pro_band_elements_mat,
            element_type, cfg)

    if "csv" in output_formats:
        print("Writing CSV files...")
        writers.write_csv_band_structure(kvec_mat, eigen_val_mat, cfg)
        writers.write_csv_occ_band(
            kvec_mat, eigen_val_mat, elec_occ_mat, cfg)
        writers.write_csv_pband_per_ion(
            kvec_mat, eigen_val_mat, pro_band_mat, cfg)
        writers.write_csv_pband_per_element(
            kvec_mat, eigen_val_mat, pro_band_elements_mat,
            element_type, cfg)

    if "matplotlib" in output_formats:
        print("Writing matplotlib plots...")
        writers.plot_band_structure(kvec_mat, eigen_val_mat, cfg)
        writers.plot_occ_band(
            kvec_mat, eigen_val_mat, elec_occ_mat, cfg)
        writers.plot_pband_per_ion(
            kvec_mat, eigen_val_mat, pro_band_mat, cfg)
        writers.plot_pband_per_element(
            kvec_mat, eigen_val_mat, pro_band_elements_mat,
            element_type, cfg)

    print("Done.")


if __name__ == "__main__":
    main()
