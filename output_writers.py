#!/usr/bin/env python3

"""Output writers for VASP band structure data.

Supports three output formats:
  - Igor Pro .itx files (preserves exact legacy format)
  - CSV files
  - Matplotlib plots (PNG/PDF/SVG)
"""

import csv
import os

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


# ---------------------------------------------------------------------------
# Igor-related helpers (moved from YK_vasp.py)
# ---------------------------------------------------------------------------

def turn_G_into_Gamma(axis_let_mat):
    """Convert 'G' labels to Igor Pro Symbol-font Gamma characters."""
    result = list(axis_let_mat)
    for i in range(len(result)):
        if result[i] == 'G':
            result[i] = " \\F'Symbol'" + result[i] + " "
        else:
            result[i] = " " + result[i] + "\\F'Symbol' "
    return result


# ---------------------------------------------------------------------------
# Igor Pro shared helpers
# ---------------------------------------------------------------------------

def _write_igor_append_traces(f, wave_names, kvec_name, starting_band):
    """Write AppendToGraph commands, batching traces in groups of 5.

    This preserves the exact trace-batching-by-5 behavior of the original.
    """
    count_band = len(wave_names)
    count_start = 0
    count_end = min(5, count_band)

    while count_band > 5:
        f.write("X AppendToGraph\t")
        for y in range(count_start, count_end):
            f.write("%s\t" % wave_names[y])
        f.write("vs\t")
        f.write("%s\t" % kvec_name)
        f.write(";DelayUpdate\n")
        count_band -= 5
        count_start += 5
        count_end += 5

    if count_start < len(wave_names):
        f.write("X AppendToGraph\t")
        for y in range(count_start, len(wave_names)):
            f.write("%s\t" % wave_names[y])
        f.write("vs\t")
        f.write("%s\t" % kvec_name)
        f.write(";DelayUpdate\n")


def _write_igor_graph_formatting(f, energy_range, cfg):
    """Write the 14-line ModifyGraph block (previously copy-pasted 4x)."""
    width = cfg.get("igor_width", 226.772)
    height = cfg.get("igor_height", 340.157)
    font_size = cfg.get("igor_font_size", 24)
    emin, emax = energy_range

    f.write("X ModifyGraph width=%.3f,height=%.3f;DelayUpdate\n" % (width, height))
    f.write("X ModifyGraph tick=2;DelayUpdate\n")
    f.write("X ModifyGraph mirror=1;DelayUpdate\n")
    f.write("X ModifyGraph nticks(left)=5;DelayUpdate\n")
    f.write("X ModifyGraph fSize=%d;DelayUpdate\n" % font_size)
    f.write("X ModifyGraph lblMargin=15;DelayUpdate\n")
    f.write("X ModifyGraph standoff=0;DelayUpdate\n")
    f.write("X ModifyGraph axisOnTop=1;DelayUpdate\n")
    f.write("X ModifyGraph axThick=2;DelayUpdate\n")
    f.write("X ModifyGraph zero(left)=8,zeroThick(left)=2;DelayUpdate\n")
    f.write("X ModifyGraph grid(bottom)=1,gridHair(bottom)=0,gridStyle(bottom)=5,gridRGB(bottom)=(0,0,0);DelayUpdate\n")
    f.write('X Label left "Energy (eV)";DelayUpdate\n')
    f.write('X Label bottom "Wavevector";DelayUpdate\n')
    f.write("X SetAxis left %g,%g\n" % (emin, emax))


# ---------------------------------------------------------------------------
# Igor Pro writers
# ---------------------------------------------------------------------------

def write_igor_band_structure(kvec, eigen_val, cfg):
    """Write band_structure.itx — plain band structure."""
    suffix = cfg["wave_suffix"]
    s_kp = cfg["starting_kpoint"] - 1
    e_kp = cfg["ending_kpoint"]
    s_band = cfg["starting_band"] - 1
    e_band = cfg["ending_band"]

    with open("band_structure.itx", "w") as f:
        f.write("IGOR\n")
        f.write("WAVES/D\t")
        f.write("Kvector_%s\t" % suffix)
        for y in range(s_band, e_band):
            f.write("E_band_%03d_%s\t" % (y + 1, suffix))
        f.write("\n")
        f.write("BEGIN\n")

        for x in range(s_kp, e_kp):
            f.write("%.10f\t" % kvec[x])
            for y in range(s_band, e_band):
                f.write("%.8f\t" % eigen_val[y][x])
            f.write("\n")

        f.write("END\n")

        # Display command
        f.write("X Display\t")
        f.write("as\t")
        f.write('"band_structure_%s"' % suffix)
        f.write(";DelayUpdate\n")

        # AppendToGraph in batches of 5
        wave_names = ["E_band_%03d_%s" % (y + 1, suffix)
                      for y in range(s_band, e_band)]
        kvec_name = "Kvector_%s" % suffix
        _write_igor_append_traces(f, wave_names, kvec_name, s_band)

        # Per-band styling
        for y in range(s_band, e_band):
            f.write("X ModifyGraph rgb(E_band_%03d_%s)=(65535,43690,0),"
                    "lsize(E_band_%03d_%s)=2;DelayUpdate\n"
                    % (y + 1, suffix, y + 1, suffix))

        _write_igor_graph_formatting(f, cfg["energy_range"], cfg)


def write_igor_occ_band(kvec, eigen_val, elec_occ, cfg):
    """Write oband_structure.itx — occupation-colored bands."""
    suffix = cfg["wave_suffix"]
    s_kp = cfg["starting_kpoint"] - 1
    e_kp = cfg["ending_kpoint"]
    s_band = cfg["starting_band"] - 1
    e_band = cfg["ending_band"]

    with open("oband_structure.itx", "w") as f:
        f.write("IGOR\n")
        f.write("WAVES/D\t")
        f.write("Kvector_occ_%s\t" % suffix)
        for y in range(s_band, e_band):
            f.write("E_band_%03d_occ_%s\t" % (y + 1, suffix))
        for y in range(s_band, e_band):
            f.write("O_band_%03d_occ_%s\t" % (y + 1, suffix))
        f.write("\n")
        f.write("BEGIN\n")

        for x in range(s_kp, e_kp):
            f.write("%.10f\t" % kvec[x])
            for y in range(s_band, e_band):
                f.write("%.8f\t" % eigen_val[y][x])
            for y in range(s_band, e_band):
                f.write("%.8f\t" % elec_occ[y][x])
            f.write("\n")

        f.write("END\n")

        f.write("X Display\t")
        f.write("as\t")
        f.write('"oband_structure_occ_%s"' % suffix)
        f.write(";DelayUpdate\n")

        wave_names = ["E_band_%03d_occ_%s" % (y + 1, suffix)
                      for y in range(s_band, e_band)]
        kvec_name = "Kvector_occ_%s" % suffix
        _write_igor_append_traces(f, wave_names, kvec_name, s_band)

        for y in range(s_band, e_band):
            f.write("X ModifyGraph mode(E_band_%03d_occ_%s)=0,lsize=2,"
                    "zColor(E_band_%03d_occ_%s)="
                    "{O_band_%03d_occ_%s,0,1,BlueGreenOrange,1}"
                    ";DelayUpdate\n"
                    % (y + 1, suffix, y + 1, suffix, y + 1, suffix))

        _write_igor_graph_formatting(f, cfg["energy_range"], cfg)


def write_igor_pband_per_ion(kvec, eigen_val, pro_band, cfg):
    """Write per-ion projected band .itx files."""
    suffix = cfg["wave_suffix"]
    s_kp = cfg["starting_kpoint"] - 1
    e_kp = cfg["ending_kpoint"]
    s_band = cfg["starting_band"] - 1
    e_band = cfg["ending_band"]
    s_ion = cfg["starting_ion"] - 1
    e_ion = cfg["ending_ion"]
    orbital_types = cfg["orbital_types"]
    proj_max = cfg.get("proj_color_max", 0.6)

    for t in range(len(orbital_types)):
        for z in range(s_ion, e_ion):
            filename = "%02dion_pband_%s.itx" % (z + 1, orbital_types[t])
            with open(filename, "w") as f:
                f.write("IGOR\n")
                f.write("WAVES/D\t")
                kvec_name = "Kvector_%02dion_%s_%s" % (
                    z + 1, orbital_types[t], suffix)
                f.write("%s\t" % kvec_name)

                e_names = []
                p_names = []
                for y in range(s_band, e_band):
                    name = "E_band_%03d_%02dion_%s_%s" % (
                        y + 1, z + 1, orbital_types[t], suffix)
                    e_names.append(name)
                    f.write("%s\t" % name)
                for y in range(s_band, e_band):
                    name = "P_band_%03d_%02dion_%s_%s" % (
                        y + 1, z + 1, orbital_types[t], suffix)
                    p_names.append(name)
                    f.write("%s\t" % name)

                f.write("\n")
                f.write("BEGIN\n")

                for x in range(s_kp, e_kp):
                    f.write("%.10f\t" % kvec[x])
                    for y in range(s_band, e_band):
                        f.write("%.8f\t" % eigen_val[y][x])
                    for y in range(s_band, e_band):
                        f.write("%.8f\t" % pro_band[t][z][y][x])
                    f.write("\n")

                f.write("END\n")

                f.write("X Display\t")
                f.write("as\t")
                f.write('"pband_%02dion_%s_%s"' % (
                    z + 1, orbital_types[t], suffix))
                f.write(";DelayUpdate\n")

                _write_igor_append_traces(f, e_names, kvec_name, s_band)

                for y in range(s_band, e_band):
                    e_n = "E_band_%03d_%02dion_%s_%s" % (
                        y + 1, z + 1, orbital_types[t], suffix)
                    p_n = "P_band_%03d_%02dion_%s_%s" % (
                        y + 1, z + 1, orbital_types[t], suffix)
                    f.write("X ModifyGraph mode(%s)=0,lsize=2,"
                            "zColor(%s)={%s,0,%g,Red,1}"
                            ";DelayUpdate\n" % (e_n, e_n, p_n, proj_max))

                _write_igor_graph_formatting(f, cfg["energy_range"], cfg)


def write_igor_pband_per_element(kvec, eigen_val, pro_band_elements,
                                 element_type, cfg):
    """Write per-element projected band .itx files."""
    suffix = cfg["wave_suffix"]
    s_kp = cfg["starting_kpoint"] - 1
    e_kp = cfg["ending_kpoint"]
    s_band = cfg["starting_band"] - 1
    e_band = cfg["ending_band"]
    orbital_types = cfg["orbital_types"]
    proj_max = cfg.get("proj_color_max", 0.6)

    for t in range(len(orbital_types)):
        for z in range(len(element_type)):
            elem = element_type[z]
            filename = "%s_pband_%s.itx" % (elem, orbital_types[t])
            with open(filename, "w") as f:
                f.write("IGOR\n")
                f.write("WAVES/D\t")
                kvec_name = "Kvector_%s_%s_%s" % (
                    elem, orbital_types[t], suffix)
                f.write("%s\t" % kvec_name)

                e_names = []
                p_names = []
                for y in range(s_band, e_band):
                    name = "E_band_%03d_%s_%s_%s" % (
                        y + 1, elem, orbital_types[t], suffix)
                    e_names.append(name)
                    f.write("%s\t" % name)
                for y in range(s_band, e_band):
                    name = "P_band_%03d_%s_%s_%s" % (
                        y + 1, elem, orbital_types[t], suffix)
                    p_names.append(name)
                    f.write("%s\t" % name)

                f.write("\n")
                f.write("BEGIN\n")

                for x in range(s_kp, e_kp):
                    f.write("%.10f\t" % kvec[x])
                    for y in range(s_band, e_band):
                        f.write("%.8f\t" % eigen_val[y][x])
                    for y in range(s_band, e_band):
                        f.write("%.8f\t" % pro_band_elements[t][z][y][x])
                    f.write("\n")

                f.write("END\n")

                f.write("X Display\t")
                f.write("as\t")
                f.write('"pband_%s_%s_%s"' % (
                    elem, orbital_types[t], suffix))
                f.write(";DelayUpdate\n")

                _write_igor_append_traces(f, e_names, kvec_name, s_band)

                for y in range(s_band, e_band):
                    e_n = "E_band_%03d_%s_%s_%s" % (
                        y + 1, elem, orbital_types[t], suffix)
                    p_n = "P_band_%03d_%s_%s_%s" % (
                        y + 1, elem, orbital_types[t], suffix)
                    f.write("X ModifyGraph mode(%s)=0,lsize=2,"
                            "zColor(%s)={%s,0,%g,Red,1}"
                            ";DelayUpdate\n" % (e_n, e_n, p_n, proj_max))

                _write_igor_graph_formatting(f, cfg["energy_range"], cfg)


def write_igor_pband_select_atoms(kvec, eigen_val, pro_band_selec, cfg):
    """Write select-atoms projected band .itx files.

    Uses "SL" prefix and BlackBody colormap, matching the original
    post_PROCAR_select_atom.py output.
    """
    suffix = cfg["wave_suffix"]
    s_kp = cfg["starting_kpoint"] - 1
    e_kp = cfg["ending_kpoint"]
    s_band = cfg["starting_band"] - 1
    e_band = cfg["ending_band"]
    orbital_types = cfg["orbital_types"]
    proj_max = cfg.get("proj_color_max", 0.6)

    for t in range(len(orbital_types)):
        orb = orbital_types[t]
        filename = "SL_%s_pband_%s.itx" % (suffix, orb)
        with open(filename, "w") as f:
            f.write("IGOR\n")
            f.write("WAVES/D\t")
            kvec_name = "Kvector_sl_%s_%s" % (orb, suffix)
            f.write("%s\t" % kvec_name)

            e_names = []
            p_names = []
            for y in range(s_band, e_band):
                name = "E_band_%03d_sl_%s_%s" % (y + 1, orb, suffix)
                e_names.append(name)
                f.write("%s\t" % name)
            for y in range(s_band, e_band):
                name = "P_band_%03d_sl_%s_%s" % (y + 1, orb, suffix)
                p_names.append(name)
                f.write("%s\t" % name)

            f.write("\n")
            f.write("BEGIN\n")

            for x in range(s_kp, e_kp):
                f.write("%.10f\t" % kvec[x])
                for y in range(s_band, e_band):
                    f.write("%.8f\t" % eigen_val[y][x])
                for y in range(s_band, e_band):
                    f.write("%.8f\t" % pro_band_selec[t][y][x])
                f.write("\n")

            f.write("END\n")

            f.write("X Display\t")
            f.write("as\t")
            f.write('"pband_sl_%s_%s"' % (orb, suffix))
            f.write(";DelayUpdate\n")

            _write_igor_append_traces(f, e_names, kvec_name, s_band)

            for y in range(s_band, e_band):
                e_n = "E_band_%03d_sl_%s_%s" % (y + 1, orb, suffix)
                p_n = "P_band_%03d_sl_%s_%s" % (y + 1, orb, suffix)
                f.write("X ModifyGraph mode(%s)=0,lsize=2,"
                        "zColor(%s)={%s,0,%g,BlackBody,1}"
                        ";DelayUpdate\n" % (e_n, e_n, p_n, proj_max))

            _write_igor_graph_formatting(f, cfg["energy_range"], cfg)


# ---------------------------------------------------------------------------
# CSV writers
# ---------------------------------------------------------------------------

def write_csv_band_structure(kvec, eigen_val, cfg):
    """Write band_structure.csv."""
    s_kp = cfg["starting_kpoint"] - 1
    e_kp = cfg["ending_kpoint"]
    s_band = cfg["starting_band"] - 1
    e_band = cfg["ending_band"]

    with open("band_structure.csv", "w", newline="") as f:
        writer = csv.writer(f)
        header = ["kvector"] + ["E_band_%03d" % (y + 1)
                                for y in range(s_band, e_band)]
        writer.writerow(header)
        for x in range(s_kp, e_kp):
            row = [kvec[x]] + [eigen_val[y][x] for y in range(s_band, e_band)]
            writer.writerow(row)


def write_csv_occ_band(kvec, eigen_val, elec_occ, cfg):
    """Write oband_structure.csv."""
    s_kp = cfg["starting_kpoint"] - 1
    e_kp = cfg["ending_kpoint"]
    s_band = cfg["starting_band"] - 1
    e_band = cfg["ending_band"]

    with open("oband_structure.csv", "w", newline="") as f:
        writer = csv.writer(f)
        header = (["kvector"]
                  + ["E_band_%03d" % (y + 1) for y in range(s_band, e_band)]
                  + ["O_band_%03d" % (y + 1) for y in range(s_band, e_band)])
        writer.writerow(header)
        for x in range(s_kp, e_kp):
            row = ([kvec[x]]
                   + [eigen_val[y][x] for y in range(s_band, e_band)]
                   + [elec_occ[y][x] for y in range(s_band, e_band)])
            writer.writerow(row)


def write_csv_pband_per_ion(kvec, eigen_val, pro_band, cfg):
    """Write per-ion projected band CSV files."""
    s_kp = cfg["starting_kpoint"] - 1
    e_kp = cfg["ending_kpoint"]
    s_band = cfg["starting_band"] - 1
    e_band = cfg["ending_band"]
    s_ion = cfg["starting_ion"] - 1
    e_ion = cfg["ending_ion"]
    orbital_types = cfg["orbital_types"]

    for t in range(len(orbital_types)):
        for z in range(s_ion, e_ion):
            filename = "%02dion_pband_%s.csv" % (z + 1, orbital_types[t])
            with open(filename, "w", newline="") as f:
                writer = csv.writer(f)
                header = (["kvector"]
                          + ["E_band_%03d" % (y + 1)
                             for y in range(s_band, e_band)]
                          + ["P_band_%03d" % (y + 1)
                             for y in range(s_band, e_band)])
                writer.writerow(header)
                for x in range(s_kp, e_kp):
                    row = ([kvec[x]]
                           + [eigen_val[y][x] for y in range(s_band, e_band)]
                           + [pro_band[t][z][y][x]
                              for y in range(s_band, e_band)])
                    writer.writerow(row)


def write_csv_pband_per_element(kvec, eigen_val, pro_band_elements,
                                element_type, cfg):
    """Write per-element projected band CSV files."""
    s_kp = cfg["starting_kpoint"] - 1
    e_kp = cfg["ending_kpoint"]
    s_band = cfg["starting_band"] - 1
    e_band = cfg["ending_band"]
    orbital_types = cfg["orbital_types"]

    for t in range(len(orbital_types)):
        for z in range(len(element_type)):
            elem = element_type[z]
            filename = "%s_pband_%s.csv" % (elem, orbital_types[t])
            with open(filename, "w", newline="") as f:
                writer = csv.writer(f)
                header = (["kvector"]
                          + ["E_band_%03d" % (y + 1)
                             for y in range(s_band, e_band)]
                          + ["P_band_%03d" % (y + 1)
                             for y in range(s_band, e_band)])
                writer.writerow(header)
                for x in range(s_kp, e_kp):
                    row = ([kvec[x]]
                           + [eigen_val[y][x] for y in range(s_band, e_band)]
                           + [pro_band_elements[t][z][y][x]
                              for y in range(s_band, e_band)])
                    writer.writerow(row)


def write_csv_pband_select_atoms(kvec, eigen_val, pro_band_selec, cfg):
    """Write select-atoms projected band CSV files."""
    s_kp = cfg["starting_kpoint"] - 1
    e_kp = cfg["ending_kpoint"]
    s_band = cfg["starting_band"] - 1
    e_band = cfg["ending_band"]
    suffix = cfg["wave_suffix"]
    orbital_types = cfg["orbital_types"]

    for t in range(len(orbital_types)):
        orb = orbital_types[t]
        filename = "SL_%s_pband_%s.csv" % (suffix, orb)
        with open(filename, "w", newline="") as f:
            writer = csv.writer(f)
            header = (["kvector"]
                      + ["E_band_%03d" % (y + 1)
                         for y in range(s_band, e_band)]
                      + ["P_band_%03d" % (y + 1)
                         for y in range(s_band, e_band)])
            writer.writerow(header)
            for x in range(s_kp, e_kp):
                row = ([kvec[x]]
                       + [eigen_val[y][x] for y in range(s_band, e_band)]
                       + [pro_band_selec[t][y][x]
                          for y in range(s_band, e_band)])
                writer.writerow(row)


# ---------------------------------------------------------------------------
# Matplotlib writers
# ---------------------------------------------------------------------------

def _check_matplotlib():
    if not HAS_MATPLOTLIB:
        raise ImportError(
            "matplotlib is required for plot output. "
            "Install with: pip install matplotlib"
        )


def _setup_band_axes(ax, energy_range):
    """Apply common band-structure axis formatting."""
    ax.set_ylabel("Energy (eV)")
    ax.set_xlabel("Wavevector")
    ax.set_ylim(energy_range)
    ax.axhline(0, color="black", linewidth=0.8, linestyle="--")


def plot_band_structure(kvec, eigen_val, cfg):
    """Plot plain band structure and save to file."""
    _check_matplotlib()
    s_kp = cfg["starting_kpoint"] - 1
    e_kp = cfg["ending_kpoint"]
    s_band = cfg["starting_band"] - 1
    e_band = cfg["ending_band"]
    fmt = cfg.get("plot_format", "png")

    k = kvec[s_kp:e_kp]

    fig, ax = plt.subplots(figsize=(4, 6))
    for y in range(s_band, e_band):
        energies = [eigen_val[y][x] for x in range(s_kp, e_kp)]
        ax.plot(k, energies, color="#FFA500", linewidth=1.0)

    _setup_band_axes(ax, cfg["energy_range"])
    fig.tight_layout()
    fig.savefig("band_structure.%s" % fmt, dpi=300)
    plt.close(fig)


def plot_occ_band(kvec, eigen_val, elec_occ, cfg):
    """Plot occupation-colored band structure."""
    _check_matplotlib()
    s_kp = cfg["starting_kpoint"] - 1
    e_kp = cfg["ending_kpoint"]
    s_band = cfg["starting_band"] - 1
    e_band = cfg["ending_band"]
    fmt = cfg.get("plot_format", "png")

    k = kvec[s_kp:e_kp]

    fig, ax = plt.subplots(figsize=(4, 6))
    for y in range(s_band, e_band):
        energies = [eigen_val[y][x] for x in range(s_kp, e_kp)]
        occ = [elec_occ[y][x] for x in range(s_kp, e_kp)]

        points = np.array([k, energies]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        # Average occupation for each segment
        seg_occ = [(occ[i] + occ[i + 1]) / 2 for i in range(len(occ) - 1)]

        lc = LineCollection(segments, cmap="RdYlBu", norm=plt.Normalize(0, 1))
        lc.set_array(np.array(seg_occ))
        lc.set_linewidth(1.5)
        ax.add_collection(lc)

    _setup_band_axes(ax, cfg["energy_range"])
    ax.set_xlim(k[0], k[-1])
    cbar = fig.colorbar(lc, ax=ax)
    cbar.set_label("Occupation")
    fig.tight_layout()
    fig.savefig("oband_structure.%s" % fmt, dpi=300)
    plt.close(fig)


def plot_pband_per_ion(kvec, eigen_val, pro_band, cfg):
    """Plot per-ion projected band structure."""
    _check_matplotlib()
    s_kp = cfg["starting_kpoint"] - 1
    e_kp = cfg["ending_kpoint"]
    s_band = cfg["starting_band"] - 1
    e_band = cfg["ending_band"]
    s_ion = cfg["starting_ion"] - 1
    e_ion = cfg["ending_ion"]
    orbital_types = cfg["orbital_types"]
    proj_max = cfg.get("proj_color_max", 0.6)
    fmt = cfg.get("plot_format", "png")

    k = kvec[s_kp:e_kp]

    for t in range(len(orbital_types)):
        for z in range(s_ion, e_ion):
            fig, ax = plt.subplots(figsize=(4, 6))
            for y in range(s_band, e_band):
                energies = [eigen_val[y][x] for x in range(s_kp, e_kp)]
                proj = [pro_band[t][z][y][x] for x in range(s_kp, e_kp)]

                points = np.array([k, energies]).T.reshape(-1, 1, 2)
                segments = np.concatenate([points[:-1], points[1:]], axis=1)
                seg_proj = [(proj[i] + proj[i + 1]) / 2
                            for i in range(len(proj) - 1)]

                lc = LineCollection(segments, cmap="Reds",
                                    norm=plt.Normalize(0, proj_max))
                lc.set_array(np.array(seg_proj))
                lc.set_linewidth(1.5)
                ax.add_collection(lc)

            _setup_band_axes(ax, cfg["energy_range"])
            ax.set_xlim(k[0], k[-1])
            cbar = fig.colorbar(lc, ax=ax)
            cbar.set_label("Projection weight")
            ax.set_title("Ion %02d — %s" % (z + 1, orbital_types[t]))
            fig.tight_layout()
            fig.savefig("%02dion_pband_%s.%s" % (z + 1, orbital_types[t], fmt),
                        dpi=300)
            plt.close(fig)


def plot_pband_per_element(kvec, eigen_val, pro_band_elements,
                           element_type, cfg):
    """Plot per-element projected band structure."""
    _check_matplotlib()
    s_kp = cfg["starting_kpoint"] - 1
    e_kp = cfg["ending_kpoint"]
    s_band = cfg["starting_band"] - 1
    e_band = cfg["ending_band"]
    orbital_types = cfg["orbital_types"]
    proj_max = cfg.get("proj_color_max", 0.6)
    fmt = cfg.get("plot_format", "png")

    k = kvec[s_kp:e_kp]

    for t in range(len(orbital_types)):
        for z in range(len(element_type)):
            elem = element_type[z]
            fig, ax = plt.subplots(figsize=(4, 6))
            for y in range(s_band, e_band):
                energies = [eigen_val[y][x] for x in range(s_kp, e_kp)]
                proj = [pro_band_elements[t][z][y][x]
                        for x in range(s_kp, e_kp)]

                points = np.array([k, energies]).T.reshape(-1, 1, 2)
                segments = np.concatenate([points[:-1], points[1:]], axis=1)
                seg_proj = [(proj[i] + proj[i + 1]) / 2
                            for i in range(len(proj) - 1)]

                lc = LineCollection(segments, cmap="Reds",
                                    norm=plt.Normalize(0, proj_max))
                lc.set_array(np.array(seg_proj))
                lc.set_linewidth(1.5)
                ax.add_collection(lc)

            _setup_band_axes(ax, cfg["energy_range"])
            ax.set_xlim(k[0], k[-1])
            cbar = fig.colorbar(lc, ax=ax)
            cbar.set_label("Projection weight")
            ax.set_title("%s — %s" % (elem, orbital_types[t]))
            fig.tight_layout()
            fig.savefig("%s_pband_%s.%s" % (elem, orbital_types[t], fmt),
                        dpi=300)
            plt.close(fig)


def plot_pband_select_atoms(kvec, eigen_val, pro_band_selec, cfg):
    """Plot select-atoms projected band structure.

    Uses 'hot' colormap to approximate Igor Pro's BlackBody.
    """
    _check_matplotlib()
    s_kp = cfg["starting_kpoint"] - 1
    e_kp = cfg["ending_kpoint"]
    s_band = cfg["starting_band"] - 1
    e_band = cfg["ending_band"]
    suffix = cfg["wave_suffix"]
    orbital_types = cfg["orbital_types"]
    proj_max = cfg.get("proj_color_max", 0.6)
    fmt = cfg.get("plot_format", "png")

    k = kvec[s_kp:e_kp]

    for t in range(len(orbital_types)):
        orb = orbital_types[t]
        fig, ax = plt.subplots(figsize=(4, 6))
        for y in range(s_band, e_band):
            energies = [eigen_val[y][x] for x in range(s_kp, e_kp)]
            proj = [pro_band_selec[t][y][x] for x in range(s_kp, e_kp)]

            points = np.array([k, energies]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            seg_proj = [(proj[i] + proj[i + 1]) / 2
                        for i in range(len(proj) - 1)]

            lc = LineCollection(segments, cmap="hot",
                                norm=plt.Normalize(0, proj_max))
            lc.set_array(np.array(seg_proj))
            lc.set_linewidth(1.5)
            ax.add_collection(lc)

        _setup_band_axes(ax, cfg["energy_range"])
        ax.set_xlim(k[0], k[-1])
        cbar = fig.colorbar(lc, ax=ax)
        cbar.set_label("Projection weight")
        ax.set_title("Select atoms — %s" % orb)
        fig.tight_layout()
        fig.savefig("SL_%s_pband_%s.%s" % (suffix, orb, fmt), dpi=300)
        plt.close(fig)
