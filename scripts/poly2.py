#!/usr/bin/env python3

"""
@file poly2.py
@author Andrea Di Antonio (github.com/diantonioandrea)
@date 2024-10-21

@copyright Copyright (c) 2024
"""

import matplotlib.pyplot as plt
import sys

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

# Black.
black: list[float] = [7 / 255, 54 / 255, 66 / 255]

# Plot.
fig, ax = plt.subplots()

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

    # Plots.

    # Lower base.
    ax.fill(x, y, facecolor=(1, 1, 1, 1), edgecolor=black, linewidth=0.25)

# Aspect.
ax.set_aspect('equal', adjustable='box')

# Output.
plt.show()