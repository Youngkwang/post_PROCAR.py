#!/usr/bin/env python3

"""Configuration module for VASP post-processing scripts.

Provides 3-layer configuration merging: DEFAULTS < YAML file < CLI args.
"""

import argparse
import os
import sys

try:
    import yaml
    HAS_YAML = True
except ImportError:
    HAS_YAML = False


DEFAULTS = {
    # Wave suffix for Igor wave names
    "wave_suffix": "default",

    # Reference energy (eV). Set to "auto" to read E-fermi from OUTCAR.
    "ref_energy": "auto",

    # K-point range (1-indexed). None means use full range.
    "starting_kpoint": 1,
    "ending_kpoint": None,

    # Band range (1-indexed). None means use full range.
    "starting_band": 1,
    "ending_band": None,

    # Ion range (1-indexed). None means use full range.
    "starting_ion": 1,
    "ending_ion": None,

    # Input file paths
    "procar_file": "PROCAR",
    "outcar_file": "OUTCAR",
    "contcar_file": "CONTCAR",

    # Output formats: any combination of "igor", "matplotlib", "csv"
    "output_formats": ["igor"],

    # Energy range for plots [min, max] in eV
    "energy_range": [-3, 3],

    # Matplotlib output format
    "plot_format": "png",

    # Igor Pro graph dimensions (points)
    "igor_width": 226.772,
    "igor_height": 340.157,

    # Igor Pro font size
    "igor_font_size": 24,

    # Orbital types for projected band output
    "orbital_types": ["1_s", "2_p", "3_d", "4_t"],

    # Projection color scale max for per-ion/element projected bands
    "proj_color_max": 0.6,
}


def build_parser():
    """Build argparse parser with all CLI flags."""
    parser = argparse.ArgumentParser(
        description="VASP post-processing: parse PROCAR/OUTCAR and generate "
                    "band structure outputs (Igor .itx, matplotlib, CSV).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
examples:
  python post_procar.py -w GaAs
  python post_procar.py -w GaAs -o igor matplotlib csv
  python post_procar.py --config my_calc.yaml
  python post_procar.py --config my_calc.yaml --ref-energy -0.5 --energy-range -5 5
""",
    )

    parser.add_argument(
        "--config", metavar="FILE",
        help="YAML configuration file (values override defaults, "
             "CLI args override config file)",
    )
    parser.add_argument(
        "-w", "--wave-suffix", metavar="NAME",
        help="Suffix for Igor wave names (default: %(default)s)",
    )
    parser.add_argument(
        "--ref-energy", type=float, metavar="EV",
        help="Reference energy in eV. Omit or set to 'auto' in YAML "
             "to auto-detect E-fermi from OUTCAR.",
    )
    parser.add_argument(
        "--starting-kpoint", type=int, metavar="N",
        help="First k-point index (1-indexed, default: 1)",
    )
    parser.add_argument(
        "--ending-kpoint", type=int, metavar="N",
        help="Last k-point index (1-indexed, default: all)",
    )
    parser.add_argument(
        "--starting-band", type=int, metavar="N",
        help="First band index (1-indexed, default: 1)",
    )
    parser.add_argument(
        "--ending-band", type=int, metavar="N",
        help="Last band index (1-indexed, default: all)",
    )
    parser.add_argument(
        "--starting-ion", type=int, metavar="N",
        help="First ion index (1-indexed, default: 1)",
    )
    parser.add_argument(
        "--ending-ion", type=int, metavar="N",
        help="Last ion index (1-indexed, default: all)",
    )
    parser.add_argument(
        "-o", "--output-formats", nargs="+",
        choices=["igor", "matplotlib", "csv"],
        help="Output formats (default: igor)",
    )
    parser.add_argument(
        "--energy-range", type=float, nargs=2, metavar=("MIN", "MAX"),
        help="Energy axis range in eV (default: -3 3)",
    )
    parser.add_argument(
        "--plot-format", choices=["png", "pdf", "svg"],
        help="Matplotlib output format (default: png)",
    )
    parser.add_argument(
        "--procar-file", metavar="FILE",
        help="Path to PROCAR file (default: PROCAR)",
    )
    parser.add_argument(
        "--outcar-file", metavar="FILE",
        help="Path to OUTCAR file (default: OUTCAR)",
    )
    parser.add_argument(
        "--contcar-file", metavar="FILE",
        help="Path to CONTCAR file (default: CONTCAR)",
    )

    return parser


def _load_yaml(filepath):
    """Load a YAML configuration file."""
    if not HAS_YAML:
        print("Error: PyYAML is required to use --config. "
              "Install with: pip install pyyaml", file=sys.stderr)
        sys.exit(1)
    with open(filepath, "r") as f:
        data = yaml.safe_load(f)
    return data if data is not None else {}


def get_config(argv=None):
    """Build configuration by merging DEFAULTS < YAML < CLI args.

    Parameters
    ----------
    argv : list or None
        Command-line arguments (default: sys.argv[1:]).

    Returns
    -------
    dict
        Merged configuration dictionary.
    """
    parser = build_parser()
    args = parser.parse_args(argv)

    # Start with defaults
    cfg = dict(DEFAULTS)

    # Layer 2: YAML file overrides
    if args.config:
        if not os.path.isfile(args.config):
            print(f"Error: config file not found: {args.config}",
                  file=sys.stderr)
            sys.exit(1)
        yaml_cfg = _load_yaml(args.config)
        # YAML keys use underscores to match DEFAULTS
        for key, value in yaml_cfg.items():
            if key in cfg:
                cfg[key] = value
            else:
                print(f"Warning: unknown config key '{key}' in {args.config}",
                      file=sys.stderr)

    # Layer 3: CLI args override (only non-None values)
    cli_map = {
        "wave_suffix": args.wave_suffix,
        "ref_energy": args.ref_energy,
        "starting_kpoint": args.starting_kpoint,
        "ending_kpoint": args.ending_kpoint,
        "starting_band": args.starting_band,
        "ending_band": args.ending_band,
        "starting_ion": args.starting_ion,
        "ending_ion": args.ending_ion,
        "output_formats": args.output_formats,
        "energy_range": args.energy_range,
        "plot_format": args.plot_format,
        "procar_file": args.procar_file,
        "outcar_file": args.outcar_file,
        "contcar_file": args.contcar_file,
    }
    for key, value in cli_map.items():
        if value is not None:
            cfg[key] = value

    return cfg
