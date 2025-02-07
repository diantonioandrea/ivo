#!/usr/bin/env python3

"""
@file poly2.py tikz wrapper.
@author Andrea Di Antonio (github.com/diantonioandrea)
@date 2025-02-07

@copyright Copyright (c) 2024
"""

import sys
import os

# Output directory.
os.makedirs("output", exist_ok=True)

# Template.
template_file = open("templates/square.tex", "r")
template: str = template_file.read()
template_file.close()

if len(sys.argv) <= 1:
    print(f"Usage: {sys.argv[0]} /path/to/file.p2.")
    sys.exit(0)

try:
    if sys.argv[1].split(".")[-1] != "p2":
        raise

    file = open(sys.argv[1], "r+")
    lines: list[str] = file.read().split("\n")
    file.close()

except FileNotFoundError:
    print("File not found.")
    sys.exit(-1)

except:
    print("Load a .p2 file.")
    sys.exit(-1)

# Cells.
cells: list[str] = []

for line in lines:
    if not line:
        continue

    if line[0] == "@":
        continue

    # Data.
    data: list[str] = line.split(" ")
    
    # Polygon.
    x: list[float] = []
    y: list[float] = []

    try:
        for index in range(0, len(data) - 1, 3):
            x.append(float(data[index]))
            y.append(float(data[index + 1]))

    except ValueError:
        continue

    # Cell.
    xy: list[str] = []

    for _x, _y in zip(x, y):
        xy.append(f"({_x},{_y})")

    # Cells.
    cells.append("\t\t{" + "--".join(xy) + "}")

# Text.
text: str = template.replace("_POLYGONS", ",\n".join(cells))

# Output filename.
filename: str = f"output/Square_{len(cells)}.tex"

# Output.
output = open(filename, "w")
output.write(text)
output.close()