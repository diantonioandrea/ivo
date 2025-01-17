#!/usr/bin/env python3

"""
@file error21.py tikz wrapper.
@author Andrea Di Antonio (github.com/diantonioandrea)
@date 2024-10-21

@copyright Copyright (c) 2024
"""

import numpy as np

import sys
import os

# Output directory.
os.makedirs("output", exist_ok=True)

# Template.
template_file = open("templates/loglog.tex", "r")
template: str = template_file.read()
template_file.close()

if len(sys.argv) <= 1:
    print(f"Usage: {sys.argv[0]} /path/to/file.e21.")
    sys.exit(0)

try:
    if sys.argv[1].split(".")[-1] != "e21":
        raise

    file = open(sys.argv[1], "r+")
    lines: list[str] = file.read().split("\n")
    file.close()

except FileNotFoundError:
    print("File not found.")
    sys.exit(-1)

except:
    print("Load a .e21 file.")
    sys.exit(-1)

# Quantities.
dofs: list[int] = []

p: list[int] = []
q: list[int] = []

h: list[float] = []
t: list[float] = []

l2l2: list[float] = []
l2T: list[float] = []
l2h1: list[float] = []
linfl2: list[float] = []

# Reading.
for line in lines:
    if not line:
        continue

    data = line.split(" ")[-1]

    if "DoFs:" in line:
        dofs.append(int(data))

    if "p:" in line:
        p.append(int(data))
    
    if "q:" in line:
        q.append(int(data))

    if "h:" in line:
        h.append(float(data))
    
    if "t:" in line:
        t.append(float(data))
    
    if "l2l2" in line:
        l2l2.append(float(data))
    
    if "l2T" in line:
        l2T.append(float(data))
    
    if "l2h1" in line:
        l2h1.append(float(data))

    if "linfl2" in line:
        linfl2.append(float(data))

# Tests.
tests: int = len(dofs)

# Rearrangement.
h = np.array(h).ravel()
t = np.array(t).ravel()

l2l2 = np.array(l2l2).ravel()
l2T = np.array(l2T).ravel()
l2h1 = np.array(l2h1).ravel()
linfl2 = np.array(linfl2).ravel()

# Output filename.
filename: str = f"output/{p[0]}_{q[0]}_" + ("par" if sum(l2h1) > 0.0 else "hyp") + "_"

if "--l2l2" in sys.argv: # l2l2 only.

    # Text.
    text: str = template

    # Comparisons.
    l2l2hc = (h / h[-1]) ** (p[0] + 1) * l2l2[-1]
    l2l2tc = (t / t[-1]) ** (q[0] + 1) * l2l2[-1]
    l2l2fc = ((h / h[-1]) ** (p[0] + 1) + (t / t[-1]) ** (q[0] + 1)) * l2l2[-1] / 2

    # Replacements.

    # Label.
    text = text.replace(" _LHP ", f"$h^{p[0]}$", 1)

    # Plots.
    text = text.replace(" _H_COMPARISON ", " ".join(str(_) for _ in zip(dofs, l2l2hc)))
    text = text.replace(" _T_COMPARISON ", " ".join(str(_) for _ in zip(dofs, l2l2tc)))
    text = text.replace(" _HT_COMPARISON ", " ".join(str(_) for _ in zip(dofs, l2l2fc)))
    text = text.replace(" _ERROR ", " ".join(str(_) for _ in zip(dofs, l2l2)))

    # Legend.
    text = text.replace(" _LHP ", f"$h^{p[0]}$", 1)
    text = text.replace(" _LTQ ", f"$\\tau^{q[0]}$")
    text = text.replace(" _LHPTQ ", f"$h^{p[0]} + \\tau^{q[0]}$")
    text = text.replace(" _LERROR ", "$\\Norm{\\Error}_{\\SpaceLp{2}(I; \\SpaceLp{2}(\\Omega))}$")

    # Output.
    output = open(filename + "l2l2.tex", "w")
    output.write(text)
    output.close()

