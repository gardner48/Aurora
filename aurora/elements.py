"""Minimal periodic-table lookup used by Aurora's core routines.

Provides :py:func:`atomic_element`, a drop-in replacement for the subset of
``omfit_classes.utils_math.atomic_element`` that Aurora uses, so that the core
simulation, source and radiation routines do not depend on ``omfit_classes``.

For each element symbol, the table holds the atomic number Z, the mass number A
and mass of the most abundant isotope, and the element name, matching the
values returned by ``omfit_classes.utils_math.atomic_element(symbol=...)``.
"""
# MIT License
#
# Copyright (c) 2021 Francesco Sciortino
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import re

__all__ = ["atomic_element"]

# symbol: (Z, A, mass [amu], name), most abundant isotope
_ELEMENTS = {
    "H": (1, 1, 1.007825, "Hydrogen"),
    "D": (1, 2, 2.014101, "Deuterium"),
    "DT": (1, 2.5, 2.5150751, "Deuterium-Tritium"),
    "T": (1, 3, 3.0160492, "Tritium"),
    "He": (2, 4, 4.002603, "Helium"),
    "Li": (3, 7, 7.016005, "Lithium"),
    "Be": (4, 9, 9.012183, "Beryllium"),
    "B": (5, 11, 11.009305, "Boron"),
    "C": (6, 12, 12.0, "Carbon"),
    "N": (7, 14, 14.003074, "Nitrogen"),
    "O": (8, 16, 15.994915, "Oxygen"),
    "F": (9, 19, 18.998403, "Fluorine"),
    "Ne": (10, 20, 19.992439, "Neon"),
    "Na": (11, 23, 22.98977, "Sodium"),
    "Mg": (12, 24, 23.985045, "Magnesium"),
    "Al": (13, 27, 26.981541, "Aluminum"),
    "Si": (14, 28, 27.976928, "Silicon"),
    "P": (15, 31, 30.973763, "Phosphorus"),
    "S": (16, 32, 31.972072, "Sulfur"),
    "Cl": (17, 35, 34.968853, "Chlorine"),
    "Ar": (18, 40, 39.962383, "Argon"),
    "K": (19, 39, 38.963708, "Potassium"),
    "Ca": (20, 40, 39.962591, "Calcium"),
    "Sc": (21, 45, 44.955914, "Scandium"),
    "Ti": (22, 48, 47.947947, "Titanium"),
    "V": (23, 51, 50.943963, "Vanadium"),
    "Cr": (24, 52, 51.94051, "Chromium"),
    "Mn": (25, 55, 54.938046, "Manganese"),
    "Fe": (26, 56, 55.934939, "Iron"),
    "Co": (27, 59, 58.933198, "Cobalt"),
    "Ni": (28, 58, 57.935347, "Nickel"),
    "Cu": (29, 63, 62.929599, "Copper"),
    "Zn": (30, 64, 63.929145, "Zinc"),
    "Ga": (31, 69, 68.925581, "Gallium"),
    "Ge": (32, 74, 73.921179, "Germanium"),
    "As": (33, 75, 74.921596, "Arsenic"),
    "Se": (34, 80, 79.916521, "Selenium"),
    "Br": (35, 79, 78.918336, "Bromine"),
    "Kr": (36, 84, 83.911506, "Krypton"),
    "Rb": (37, 85, 84.9118, "Rubidium"),
    "Sr": (38, 88, 87.905625, "Strontium"),
    "Y": (39, 89, 88.905856, "Yttrium"),
    "Zr": (40, 90, 89.904708, "Zirconium"),
    "Nb": (41, 93, 92.906378, "Niobium"),
    "Mo": (42, 98, 97.905405, "Molybdenum"),
    "Ru": (44, 102, 101.904348, "Ruthenium"),
    "Rh": (45, 103, 102.905503, "Rhodium"),
    "Pd": (46, 106, 105.903475, "Palladium"),
    "Ag": (47, 107, 106.905095, "Silver"),
    "Cd": (48, 114, 113.903361, "Cadmium"),
    "In": (49, 115, 114.903875, "Indium"),
    "Sn": (50, 120, 119.902199, "Tin"),
    "Sb": (51, 121, 120.903824, "Antimony"),
    "Te": (52, 130, 129.906229, "Tellurium"),
    "I": (53, 127, 126.904477, "Iodine"),
    "Xe": (54, 132, 131.904148, "Xenon"),
    "Cs": (55, 133, 132.905433, "Cesium"),
    "Ba": (56, 138, 137.905236, "Barium"),
    "La": (57, 139, 138.906355, "Lanthanum"),
    "Ce": (58, 140, 139.905442, "Cerium"),
    "Pr": (59, 141, 140.907657, "Praseodymium"),
    "Nd": (60, 142, 141.907731, "Neodymium"),
    "Sm": (62, 152, 151.919741, "Samarium"),
    "Eu": (63, 153, 152.921243, "Europium"),
    "Gd": (64, 158, 157.924111, "Gadolinium"),
    "Tb": (65, 159, 158.92535, "Terbium"),
    "Dy": (66, 164, 163.929183, "Dysprosium"),
    "Ho": (67, 165, 164.930332, "Holmium"),
    "Er": (68, 166, 165.930305, "Erbium"),
    "Tm": (69, 169, 168.934225, "Thulium"),
    "Yb": (70, 174, 173.938873, "Ytterbium"),
    "Lu": (71, 175, 174.940785, "Lutetium"),
    "Hf": (72, 180, 179.946561, "Hafnium"),
    "Ta": (73, 181, 180.948014, "Tantalum"),
    "W": (74, 184, 183.950953, "Tungsten"),
    "Re": (75, 187, 186.955765, "Rhenium"),
    "Os": (76, 192, 191.961487, "Osmium"),
    "Ir": (77, 193, 192.962942, "Iridium"),
    "Pt": (78, 195, 194.964785, "Platinum"),
    "Au": (79, 197, 196.96656, "Gold"),
    "Hg": (80, 202, 201.970632, "Mercury"),
    "Tl": (81, 205, 204.97441, "Thallium"),
    "Pb": (82, 208, 207.976641, "Lead"),
    "Bi": (83, 209, 208.980388, "Bismuth"),
    "Th": (90, 232, 232.038054, "Thorium"),
    "U": (92, 238, 238.050786, "Uranium"),
}


def atomic_element(symbol):
    """Return Z and mass number information for an element symbol.

    Parameters
    ----------
    symbol : str
        Atomic symbol, e.g. ``"Ar"`` or ``"D"``. A mass number may follow the
        symbol (``"Ar40"``) or precede it, with the ion charge after it
        (``"2H1"``), as accepted by ``omfit_classes.utils_math.atomic_element``.

    Returns
    -------
    dict
        Dictionary with a single entry keyed by ``"<symbol>(<A>)"``, whose value
        has keys ``symbol``, ``symbol_A``, ``name``, ``Z``, ``Z_ion``, ``A`` and
        ``mass``. If a mass number is given in ``symbol``, ``A`` is set to it and
        ``mass`` is that of the most abundant isotope.
    """
    match = re.match(r"(\d*)([a-zA-Z]+)(\d*)$", str(symbol))
    if not match:
        raise ValueError(
            f"Wrong form of symbol, expected H, H2, or 2H1 or similar, but got {symbol}"
        )
    pre, sym, post = match.groups()
    if sym not in _ELEMENTS:
        raise ValueError(f"No atomic element satisfies symbol={sym!r}")
    Z, A, mass, name = _ELEMENTS[sym]
    Z_ion = Z
    if pre == "" and post != "":
        A = int(post)
    elif pre != "" and post != "":
        A = int(pre)
        Z_ion = int(post)
    return {
        f"{sym}({A})": {
            "symbol": sym,
            "symbol_A": f"{sym}({A})",
            "name": name,
            "Z": Z,
            "Z_ion": Z_ion,
            "A": A,
            "mass": mass,
        }
    }
