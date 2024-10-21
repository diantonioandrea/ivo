#!/usr/bin/env python3

"""
@file visual21.py
@author Andrea Di Antonio (github.com/diantonioandrea)
@date 2024-10-21

@copyright Copyright (c) 2024
"""

import matplotlib.pyplot as plt
import sys

if len(sys.argv) <= 1:
    print(f"Usage: {sys.argv[0]} /path/to/file.s21.")
    sys.exit(0)

try:
    if sys.argv[1].split(".")[-1] != "s21":
        raise

    file = open(sys.argv[1], "r+")
    lines: list[str] = file.read().split("\n")
    file.close()

except FileNotFoundError:
    print("File not found.")
    sys.exit(-1)

except:
    print("Load a .s21 file.")
    sys.exit(-1)

# Axis.
ax = plt.figure().add_subplot(projection='3d')

# Data.
x: list[float] = []
y: list[float] = []
t: list[float] = []
uh: list[float] = []

for line in lines: # Draws every polyhedron. Expensive.
    if not line:
        continue

    # Data.
    data: list[str] = line.split(",")

    x.append(float(data[0]))
    y.append(float(data[1]))
    t.append(float(data[2]))
    uh.append(float(data[3]))

# Scatter plot.
scatter = ax.scatter(x, y, t, c=uh)

# Colorbar.
plt.colorbar(scatter)

# No grid.
ax.grid(False)

# No background.
ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

plt.show()