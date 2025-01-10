#!/usr/bin/env python3

"""
@file error21.py
@author Andrea Di Antonio (github.com/diantonioandrea)
@date 2024-10-21

@copyright Copyright (c) 2024
"""

import matplotlib.pyplot as plt
import numpy as np

import sys

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

# Black.
black: list[float] = [7 / 255, 54 / 255, 66 / 255]
red: list[float] = [220 / 255, 50 / 255, 47 / 255]
green: list[float] = [133 / 255, 153 / 255, 0]

# Quantities.
dofs: list[int] = []

p: list[int] = []
q: list[int] = []

h: list[float] = []
t: list[float] = []

l2l2: list[float] = []
l2T: list[float] = []
l2h1: list[float] = []

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

# Tests.
tests: int = len(dofs)

# Rearrangement.
h = np.array(h).ravel()
t = np.array(t).ravel()

l2l2 = np.array(l2l2).ravel()
l2T = np.array(l2T).ravel()
l2h1 = np.array(l2h1).ravel()

# Plots.

# Plot flag.
pflag: bool = False

if "--l2l2" in sys.argv: # l2l2 only.

    # Figure and axes.
    fig_l2l2 = plt.figure(0)
    ax_l2l2 = fig_l2l2.gca()

    # Plot flag.
    pflag = True

    fig_l2l2.suptitle("$L^2(L^2)$ error vs. $DoFs$ on $p = " + str(p[0]) + "$ and $q = " + str(q[0]) + "$")

    # Comparisons.
    l2l2hc = (h / h[-1]) ** (p[0] + 1) * l2l2[-1]
    l2l2tc = (t / t[-1]) ** (q[0] + 1) * l2l2[-1]
    l2l2fc = ((h / h[-1]) ** (p[0] + 1) + (t / t[-1]) ** (q[0] + 1)) * l2l2[-1] / 2

    # L2(L2).
    ax_l2l2.plot(dofs, l2l2, color=black, marker="*", linewidth=1, label="$L^2(L^2)$ error")

    if p[0] == q[0]:
        ax_l2l2.plot(dofs, l2l2fc, color=green, linestyle="-", linewidth=0.75, label="$h^{" + str(p[0] + 1) + "} + \\tau^{" + str(q[0] + 1) + "}$")
    
    ax_l2l2.plot(dofs, l2l2hc, color=red, linestyle="-.", linewidth=0.5, label="$h^{" + str(p[0] + 1) + "}$")
    ax_l2l2.plot(dofs, l2l2tc, color=red, linestyle="--", linewidth=0.5, label="$\\tau^{" + str(q[0] + 1) + "}$")

    # Scale.
    ax_l2l2.set_xscale("log")
    ax_l2l2.set_yscale("log")

    # Label.
    ax_l2l2.set_xlabel("$DoFs$")

    # Legend.
    ax_l2l2.legend(loc="best")

if "--l2T" in sys.argv: # l2T.

    # Figure and axes.
    fig_l2T = plt.figure(1)
    ax_l2T = fig_l2T.gca()

    # Plot flag.
    pflag = True

    fig_l2T.suptitle("$L^2(T)$ error vs. $DoFs$ on $p = " + str(p[0]) + "$ and $q = " + str(q[0]) + "$")

    # Comparisons.
    l2Ttc = (t / t[-1]) ** (q[0] + 1) * l2T[-1]
    l2Thc = (h / h[-1]) ** (p[0] + 1) * l2T[-1]
    l2Tfc = ((h / h[-1]) ** (p[0] + 1) + (t / t[-1]) ** (q[0] + 1)) * l2T[-1] / 2

    # L2(T).
    ax_l2T.plot(dofs, l2T, color=black, marker="*", linewidth=1, label="$L^2(T)$ error")

    if p[0] == q[0]:
        ax_l2T.plot(dofs, l2Tfc, color=green, linestyle="-", linewidth=0.75, label="$h^{" + str(p[0] + 1) + "} + \\tau^{" + str(q[0] + 1) + "}$")

    ax_l2T.plot(dofs, l2Thc, color=red, linestyle="-.", linewidth=0.5, label="$h^{" + str(p[0] + 1) + "}$")
    ax_l2T.plot(dofs, l2Ttc, color=red, linestyle="--", linewidth=0.5, label="$\\tau^{" + str(q[0] + 1) + "}$")

    # Scale.
    ax_l2T.set_xscale("log")
    ax_l2T.set_yscale("log")

    # Label.
    ax_l2T.set_xlabel("$DoFs$")

    # Legend.
    ax_l2T.legend(loc="best")
    
if "--l2h1" in sys.argv: # l2h1.

    # Figure and axes.
    fig_l2h1 = plt.figure(2)
    ax_l2h1 = fig_l2h1.gca()

    # Plot flag.
    pflag = True

    fig_l2h1.suptitle("$\\sqrt{\\varepsilon}L^2(H^1)$ error vs. $DoFs$ on $p = " + str(p[0]) + "$ and $q = " + str(q[0]) + "$")

    l2h1hc = (h / h[-1]) ** p[0] * l2h1[-1]
    l2h1tc = (t / t[-1]) ** q[0] * l2h1[-1]
    l2h1fc = ((h / h[-1]) ** p[0] + (t / t[-1]) ** q[0]) * l2h1[-1] / 2

    # L2(H1).
    ax_l2h1.plot(dofs, l2h1, color=black, marker="*", linewidth=1, label="$\\sqrt{\\varepsilon}L^2(H^1)$ error")

    if p[0] == q[0]:
        ax_l2h1.plot(dofs, l2h1fc, color=green, linestyle="-", linewidth=0.75, label="$h^{" + str(p[0]) + "} + \\tau^{" + str(q[0]) + "}$")
    
    ax_l2h1.plot(dofs, l2h1hc, color=red, linestyle="-.", linewidth=0.5, label="$h^{" + str(p[0]) + "}$")
    ax_l2h1.plot(dofs, l2h1tc, color=red, linestyle="--", linewidth=0.5, label="$\\tau^{" + str(q[0]) + "}$")

    # Scale.
    ax_l2h1.set_xscale("log")
    ax_l2h1.set_yscale("log")

    # Label.
    ax_l2h1.set_xlabel("$DoFs$")

    # Legend.
    ax_l2h1.legend(loc="best")

# Output.
if pflag:
    try:
        plt.show()

    except KeyboardInterrupt:
        sys.exit(0)