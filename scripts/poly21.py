#!/usr/bin/env python3

"""
@file poly21.py
@author Andrea Di Antonio (github.com/diantonioandrea)
@date 2024-10-21

@copyright Copyright (c) 2024
"""

import matplotlib.pyplot as plt
import sys

if len(sys.argv) <= 1:
    print(f"Usage: {sys.argv[0]} /path/to/file.p21.")
    sys.exit(0)

try:
    if sys.argv[1].split(".")[-1] != "p21":
        raise

    file = open(sys.argv[1], "r+")
    lines: list[str] = file.read().split("\n")
    file.close()

except FileNotFoundError:
    print("File not found.")
    sys.exit(-1)

except:
    print("Load a .p21 file.")
    sys.exit(-1)

# Black.
black: list[float] = [7 / 255, 54 / 255, 66 / 255]

# Axis.
ax = plt.figure().add_subplot(projection='3d')

for line in lines: # Draws every polyhedron. Expensive.
    if not line:
        continue

    # Data.
    data: list[str] = line.split(",")
    
    # Lower polygon.
    lower_x: list[float] = []
    lower_y: list[float] = []
    lower_t: list[float] = []

    try:
        # Height.
        height: float = float(data[-1])

        for index in range(0, len(data) - 1, 3):
            lower_x.append(float(data[index]))
            lower_y.append(float(data[index + 1]))
            lower_t.append(float(data[index + 2]))

    except ValueError:
        continue

    # Point copy.
    lower_x.append(lower_x[0])
    lower_y.append(lower_y[0])
    lower_t.append(lower_t[0])

    # Upper polygon.
    upper_x: list[float] = lower_x
    upper_y: list[float] = lower_y
    upper_t: list[float] = [z + height for z in lower_t]

    # Plots.

    # Lower base.
    ax.plot(lower_x, lower_y, lower_t, color=black, linewidth=0.5)

    # Upper base.
    ax.plot(upper_x, upper_y, upper_t, color=black, linewidth=0.5)

    # Heights.
    for index in range(len(lower_x)):
        x: list[float] = [lower_x[index], upper_x[index]]
        y: list[float] = [lower_y[index], upper_y[index]]
        z: list[float] = [lower_t[index], upper_t[index]]

        ax.plot(x, y, z, color=black, linewidth=0.5)

# No grid.
ax.grid(False)

# No background.
ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

plt.show()