#!/usr/bin/env python3
"""Compare geometric summaries of two AMReX ASCII STL surface files."""

import argparse
import math
from pathlib import Path


def _cross(a, b, c):
    ab = tuple(b[d] - a[d] for d in range(3))
    ac = tuple(c[d] - a[d] for d in range(3))
    return (
        ab[1] * ac[2] - ab[2] * ac[1],
        ab[2] * ac[0] - ab[0] * ac[2],
        ab[0] * ac[1] - ab[1] * ac[0],
    )


def read_ascii_stl(path):
    triangles = []
    vertices = []
    solid_records = 0
    endsolid_records = 0
    facet_records = 0

    with path.open(encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, 1):
            fields = line.split()
            if not fields:
                continue
            if fields[0] == "solid":
                solid_records += 1
            elif fields[0] == "endsolid":
                endsolid_records += 1
            elif fields[:2] == ["facet", "normal"]:
                facet_records += 1
            elif fields[0] == "vertex":
                if len(fields) != 4:
                    raise ValueError(f"{path}:{line_number}: invalid vertex record")
                vertex = tuple(float(value) for value in fields[1:])
                if not all(math.isfinite(value) for value in vertex):
                    raise ValueError(f"{path}:{line_number}: non-finite vertex")
                vertices.append(vertex)
                if len(vertices) == 3:
                    triangles.append(tuple(vertices))
                    vertices.clear()

    if (solid_records != 1 or endsolid_records != 1 or vertices
            or facet_records != len(triangles)):
        raise ValueError(f"{path}: incomplete ASCII STL solid or facet records")
    if not triangles:
        raise ValueError(f"{path}: STL contains no facets")
    return triangles


def summarize(triangles, path):
    bounds_lo = [math.inf] * 3
    bounds_hi = [-math.inf] * 3
    area = 0.0
    first_moment = [0.0] * 3
    second_moment = [0.0] * 3
    orientation = [0.0] * 3

    for index, triangle in enumerate(triangles):
        for vertex in triangle:
            for d in range(3):
                bounds_lo[d] = min(bounds_lo[d], vertex[d])
                bounds_hi[d] = max(bounds_hi[d], vertex[d])

        cross = _cross(*triangle)
        twice_area = math.sqrt(sum(value * value for value in cross))
        if not twice_area > 0.0:
            raise ValueError(f"{path}: degenerate facet {index}")
        facet_area = 0.5 * twice_area
        centroid = tuple(sum(vertex[d] for vertex in triangle) / 3.0
                         for d in range(3))
        area += facet_area
        for d in range(3):
            first_moment[d] += facet_area * centroid[d]
            second_moment[d] += facet_area * centroid[d] * centroid[d]
            unit_normal = cross[d] / twice_area
            orientation[d] += facet_area * unit_normal * unit_normal

    centroid = tuple(value / area for value in first_moment)
    spread = tuple(math.sqrt(max(0.0, second_moment[d] / area
                                 - centroid[d] * centroid[d]))
                   for d in range(3))
    orientation = tuple(value / area for value in orientation)
    return {
        "facets": len(triangles),
        "area": area,
        "bounds_lo": tuple(bounds_lo),
        "bounds_hi": tuple(bounds_hi),
        "centroid": centroid,
        "spread": spread,
        "orientation": orientation,
    }


def compare(reference, candidate, relative_tolerance, absolute_tolerance):
    failures = []

    def relative_error(name, lhs, rhs):
        scale = max(abs(lhs), abs(rhs))
        error = 0.0 if scale == 0.0 else abs(lhs - rhs) / scale
        if error > relative_tolerance:
            failures.append(
                f"{name} relative error {error:.6g} exceeds "
                f"{relative_tolerance:.6g}")

    relative_error("facet count", reference["facets"], candidate["facets"])
    relative_error("surface area", reference["area"], candidate["area"])

    for name in ("bounds_lo", "bounds_hi", "centroid"):
        for d, (lhs, rhs) in enumerate(zip(reference[name], candidate[name])):
            error = abs(lhs - rhs)
            if error > absolute_tolerance:
                failures.append(
                    f"{name}[{d}] absolute error {error:.6g} exceeds "
                    f"{absolute_tolerance:.6g}")

    for name in ("spread", "orientation"):
        for d, (lhs, rhs) in enumerate(zip(reference[name], candidate[name])):
            relative_error(f"{name}[{d}]", lhs, rhs)

    return failures


def main():
    parser = argparse.ArgumentParser(
        description="Compare the shape of legacy and marching-cubes ASCII STL output")
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--relative-tolerance", type=float, default=0.02)
    parser.add_argument("--absolute-tolerance", type=float, default=0.005)
    args = parser.parse_args()

    if args.relative_tolerance < 0.0 or args.absolute_tolerance < 0.0:
        parser.error("tolerances must be non-negative")

    reference = summarize(read_ascii_stl(args.reference), args.reference)
    candidate = summarize(read_ascii_stl(args.candidate), args.candidate)
    failures = compare(reference, candidate, args.relative_tolerance,
                       args.absolute_tolerance)

    print(f"reference: facets={reference['facets']} area={reference['area']:.12g}")
    print(f"candidate: facets={candidate['facets']} area={candidate['area']:.12g}")
    if failures:
        for failure in failures:
            print(f"error: {failure}")
        raise SystemExit(1)
    print("legacy and marching-cubes STL surfaces are similar")


if __name__ == "__main__":
    main()
