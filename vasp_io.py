#!/usr/bin/env python3

"""VASP file parsing utilities.

Modernized from YK_vasp.py (2015-2017).
Parses CONTCAR, OUTCAR, and related VASP output files.
"""

import math
import numpy as np
from operator import itemgetter


# ---------------------------------------------------------------------------
# General file utilities
# ---------------------------------------------------------------------------

def get_keyword_from(keyword, filepath):
    """Return the first line containing *keyword* in *filepath*, or None."""
    with open(filepath, "r") as f:
        for line in f:
            if keyword in line:
                return line
    return None


def get_keywords_from(keyword, filepath):
    """Return all lines containing *keyword* in *filepath*."""
    with open(filepath, "r") as f:
        lines = []
        for line in f:
            if keyword in line:
                lines.append(line)
        return lines


def get_line_from(line_number, filepath):
    """Return line *line_number* (1-indexed) from *filepath*."""
    with open(filepath, "r") as f:
        target_lines = f.readlines()
        return target_lines[line_number - 1]


# ---------------------------------------------------------------------------
# CONTCAR related
# ---------------------------------------------------------------------------

def get_elements_type(target="CONTCAR"):
    """Return list of element symbols from CONTCAR line 6."""
    with open(target, "r") as f:
        lines = f.readlines()
        return lines[5].split()


def get_elements_num(target="CONTCAR"):
    """Return list of element counts (as strings) from CONTCAR line 7."""
    with open(target, "r") as f:
        lines = f.readlines()
        return lines[6].split()


def get_elements_info(target="CONTCAR"):
    """Return dict mapping element symbol → count string."""
    elements_info = {}
    types = get_elements_type(target)
    nums = get_elements_num(target)
    for i in range(len(types)):
        elements_info[types[i]] = nums[i]
    return elements_info


def get_cell_vectors(target="CONTCAR"):
    """Return dict with keys 'a', 'b', 'c' — each a list of 3 floats (Å)."""
    with open(target, "r") as f:
        lines = f.readlines()

    scale_factor = float(lines[1].split()[0])

    cell_vector_a = [float(x) * scale_factor for x in lines[2].split()]
    cell_vector_b = [float(x) * scale_factor for x in lines[3].split()]
    cell_vector_c = [float(x) * scale_factor for x in lines[4].split()]

    return {'a': cell_vector_a, 'b': cell_vector_b, 'c': cell_vector_c}


def get_lattice_constants(target="CONTCAR"):
    """Return dict with keys 'a0', 'b0', 'c0' — lattice constants in Å.

    Delegates to get_cell_vectors() instead of re-parsing the file.
    """
    cell_vectors = get_cell_vectors(target)

    def _norm(v):
        return math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)

    return {
        'a0': _norm(cell_vectors['a']),
        'b0': _norm(cell_vectors['b']),
        'c0': _norm(cell_vectors['c']),
    }


def get_atoms_position(target="CONTCAR", position_type="default"):
    """Return dict mapping element symbol → sorted list of atom positions.

    Parameters
    ----------
    target : str
        Path to CONTCAR file.
    position_type : str
        "default" for raw coordinates, "c" for Cartesian conversion.
    """
    with open(target, "r") as f:
        lines = f.readlines()

    elements_type = get_elements_type(target)
    elements_num = get_elements_num(target)
    cell_vectors = get_cell_vectors(target)

    if lines[7].split()[0][0] in ("S", "s"):
        first_atom_line = 9
    else:
        first_atom_line = 8

    atoms_position = {}
    for i in range(len(elements_type)):
        atoms_position[elements_type[i]] = [
            [0.0, 0.0, 0.0] for _ in range(int(elements_num[i]))
        ]

    element_position_line = first_atom_line

    if position_type == "default":
        for i in range(len(elements_type)):
            for j in range(element_position_line,
                           element_position_line + int(elements_num[i])):
                atoms_position[elements_type[i]][j - element_position_line] = [
                    float(lines[j].split()[k]) for k in range(3)
                ]
            element_position_line += int(elements_num[i])

    elif position_type == "c":
        coord_line = lines[first_atom_line - 1].split()[0][0]
        if coord_line in ("D", "d"):
            for i in range(len(elements_type)):
                for j in range(element_position_line,
                               element_position_line + int(elements_num[i])):
                    parts = lines[j].split()
                    pos = (np.asarray(cell_vectors["a"]) * float(parts[0])
                           + np.asarray(cell_vectors["b"]) * float(parts[1])
                           + np.asarray(cell_vectors["c"]) * float(parts[2]))
                    atoms_position[elements_type[i]][
                        j - element_position_line
                    ] = pos.tolist()
                element_position_line += int(elements_num[i])

    for elem in atoms_position:
        atoms_position[elem] = sorted(
            atoms_position[elem], key=itemgetter(2, 0, 1)
        )

    return atoms_position


# ---------------------------------------------------------------------------
# OUTCAR related
# ---------------------------------------------------------------------------

def get_VBM_band_num(target="OUTCAR"):
    """Return the VBM band number (NELEC / 2)."""
    line = get_keyword_from("NELEC", target)
    return int(float(line.split()[2]) / 2)


def get_CBM_band_num(target="OUTCAR"):
    """Return the CBM band number (VBM + 1)."""
    return get_VBM_band_num(target) + 1
