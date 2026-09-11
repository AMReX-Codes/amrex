#!/usr/bin/env python3
"""
Generate a simple ASCII STL box for EB_WMLES testing.

Usage: python make_stl.py [output_file]

The box is centered at the origin with dimensions Lx x Ly x Lz.
Adjust the parameters below to match your desired geometry.
"""

import sys, os

# ---- Box parameters (edit as needed) ----
Lx = 1.0    # length in x (streamwise)
Ly = 0.5    # length in y (spanwise / normal)
Lz = 0.5    # length in z (vertical)

outfile = sys.argv[1] if len(sys.argv) > 1 else "box.stl"

hx = Lx / 2.0
hy = Ly / 2.0
hz = Lz / 2.0

# 8 vertices of the box
v = [
    (-hx, -hy, -hz),  # 0
    ( hx, -hy, -hz),  # 1
    ( hx,  hy, -hz),  # 2
    (-hx,  hy, -hz),  # 3
    (-hx, -hy,  hz),  # 4
    ( hx, -hy,  hz),  # 5
    ( hx,  hy,  hz),  # 6
    (-hx,  hy,  hz),  # 7
]

# 12 triangles (2 per face), outward normals
# fmt: (normal, v0, v1, v2)
faces = [
    # -X face
    (( -1, 0, 0), 0, 4, 7),  (( -1, 0, 0), 0, 7, 3),
    # +X face
    ((  1, 0, 0), 1, 2, 6),  ((  1, 0, 0), 1, 6, 5),
    # -Y face
    (( 0, -1, 0), 0, 1, 5),  (( 0, -1, 0), 0, 5, 4),
    # +Y face
    (( 0,  1, 0), 2, 3, 7),  (( 0,  1, 0), 2, 7, 6),
    # -Z face
    (( 0, 0, -1), 0, 3, 2),  (( 0, 0, -1), 0, 2, 1),
    # +Z face
    (( 0, 0,  1), 4, 5, 6),  (( 0, 0,  1), 4, 6, 7),
]

def stl_vertex(vi):
    return "      vertex  {:.6f}  {:.6f}  {:.6f}".format(v[vi][0], v[vi][1], v[vi][2])

with open(outfile, "w") as f:
    f.write("solid box\n")
    for (nx, ny, nz), v0, v1, v2 in faces:
        f.write("  facet normal  {:.1f}  {:.1f}  {:.1f}\n".format(nx, ny, nz))
        f.write("    outer loop\n")
        f.write(stl_vertex(v0) + "\n")
        f.write(stl_vertex(v1) + "\n")
        f.write(stl_vertex(v2) + "\n")
        f.write("    endloop\n")
        f.write("  endfacet\n")
    f.write("endsolid box\n")

print(f"Wrote {outfile}: box {Lx} x {Ly} x {Lz}, 12 triangles")
