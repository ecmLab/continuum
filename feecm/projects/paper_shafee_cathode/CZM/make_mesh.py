#!/usr/bin/env python3
"""
make_mesh.py
============

Generate the 2-D mesh used by the cohesive-zone version of
AG_Contact_Loss.i.

Geometry (matching the original input_mesh_Shafee.msh):

      (0,L) +------------------------------+ (L,L)
            |                              |
            |                              |
            |          block_LPS           |
            |          (L-shape)           |
            |                              |
       (0,R)*\                             |
            | \. arc (block_NMC_right /    |
            |   \    block_LPS_left)       |
            |    \.                        |
            |      \.                      |
            |  block \.                    |
            |  _NMC    \.                  |
            |           \.                 |
            |             \.               |
            +-------------*----------------+
        (0,0)            (R,0)           (L,0)

    L = 5.0 mm                       (full domain side length)
    R = 4.460303300624443 mm         (NMC quarter-circle radius)

The mesh is *conformal* across the arc (shared nodes).  The MOOSE input
runs ``BreakMeshByBlockGenerator`` on it, which duplicates the interface
nodes and creates these sidesets:

    interface             # union of both sides  ->  used by the CZM action
    block_NMC_block_LPS   # NMC-side faces       ->  diagnostics
    block_LPS_block_NMC   # LPS-side faces       ->  used by NodalValueSampler

Requirements
------------
    pip install gmsh                # the Python wrapper around libgmsh

Usage
-----
    python make_mesh.py             # writes input_mesh_czm.msh
    python make_mesh.py --gui       # also open the gmsh GUI for inspection
"""

import sys
import gmsh


# ---------------------------------------------------------------------------
# Geometry / discretization parameters (taken from input_mesh_Shafee.msh)
# ---------------------------------------------------------------------------
L    = 5.0e-3                  # square domain side length [m]
R    = 4.460303300624443e-3    # NMC quarter-circle radius [m]
H_EL = 3.2e-5                  # target element edge length [m]
                               #   ~2_500 elements at 1.0e-4
                               #   ~10_000 elements at 5.0e-5
                               #   ~24_000 elements at 3.2e-5  (matches original)

# Physical-group IDs (kept consistent with the original mesh).  These are
# arbitrary in MOOSE - it matches by name - but using the same IDs makes
# diff-ing the .msh files easier.
PG = {
    "block_NMC_right": 1,   # arc (kept for backward-compat with the old .i)
    "block_LPS_left":  2,   # same arc (kept for backward-compat)
    "block_bottom":    3,
    "block_right":     4,
    "block_top":       5,
    "block_left":      6,
    "block_NMC":       7,
    "block_LPS":       8,
}

OUT_FILE = "input_mesh_czm.msh"


def build_mesh(filename: str = OUT_FILE, gui: bool = False) -> None:
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("nmc_lps_czm")

    # -----------------------------------------------------------------
    # Points
    # -----------------------------------------------------------------
    p_o   = gmsh.model.geo.addPoint(0.0, 0.0, 0.0, H_EL)   # origin / arc center
    p_Rx  = gmsh.model.geo.addPoint(R,   0.0, 0.0, H_EL)   # arc end on x-axis
    p_Ry  = gmsh.model.geo.addPoint(0.0, R,   0.0, H_EL)   # arc end on y-axis
    p_BR  = gmsh.model.geo.addPoint(L,   0.0, 0.0, H_EL)   # bottom-right corner
    p_TR  = gmsh.model.geo.addPoint(L,   L,   0.0, H_EL)   # top-right corner
    p_TL  = gmsh.model.geo.addPoint(0.0, L,   0.0, H_EL)   # top-left corner

    # -----------------------------------------------------------------
    # Curves
    #
    # NMC quarter-circle (3 boundary curves):
    #     l_nmc_bot   - along x-axis, origin -> (R, 0)
    #     l_arc       - quarter-circle arc, (R, 0) -> (0, R)
    #     l_nmc_left  - along y-axis, (0, R) -> origin
    #
    # LPS L-shape (5 boundary curves; the arc is shared):
    #     l_lps_bot   - bottom of LPS, (R, 0) -> (L, 0)
    #     l_right     - right side,    (L, 0) -> (L, L)
    #     l_top       - top side,      (L, L) -> (0, L)
    #     l_lps_left  - upper-left,    (0, L) -> (0, R)
    #     -l_arc      - shared arc, traversed in reverse
    # -----------------------------------------------------------------
    l_nmc_bot  = gmsh.model.geo.addLine(p_o,  p_Rx)
    l_arc      = gmsh.model.geo.addCircleArc(p_Rx, p_o, p_Ry)   # start, center, end
    l_nmc_left = gmsh.model.geo.addLine(p_Ry, p_o)

    l_lps_bot  = gmsh.model.geo.addLine(p_Rx, p_BR)
    l_right    = gmsh.model.geo.addLine(p_BR, p_TR)
    l_top      = gmsh.model.geo.addLine(p_TR, p_TL)
    l_lps_left = gmsh.model.geo.addLine(p_TL, p_Ry)

    # -----------------------------------------------------------------
    # Surfaces
    # -----------------------------------------------------------------
    cl_nmc = gmsh.model.geo.addCurveLoop([l_nmc_bot, l_arc, l_nmc_left])
    s_nmc  = gmsh.model.geo.addPlaneSurface([cl_nmc])

    cl_lps = gmsh.model.geo.addCurveLoop(
        [l_lps_bot, l_right, l_top, l_lps_left, -l_arc]
    )
    s_lps  = gmsh.model.geo.addPlaneSurface([cl_lps])

    gmsh.model.geo.synchronize()

    # -----------------------------------------------------------------
    # Physical groups
    #
    # block_NMC_right and block_LPS_left both point at the same arc so
    # the original AG_Contact_Loss.i (penalty-contact) can still read
    # this mesh.  For the CZM run we don't reference them - the
    # 'interface', 'block_NMC_block_LPS', and 'block_LPS_block_NMC'
    # sidesets are produced inside MOOSE by BreakMeshByBlockGenerator.
    # -----------------------------------------------------------------
    gmsh.model.addPhysicalGroup(1, [l_arc],                  PG["block_NMC_right"], "block_NMC_right")
    # gmsh.model.addPhysicalGroup(1, [l_arc],                  PG["block_LPS_left"],  "block_LPS_left")
    gmsh.model.addPhysicalGroup(1, [l_nmc_bot, l_lps_bot],   PG["block_bottom"],    "block_bottom")
    gmsh.model.addPhysicalGroup(1, [l_right],                PG["block_right"],     "block_right")
    gmsh.model.addPhysicalGroup(1, [l_top],                  PG["block_top"],       "block_top")
    gmsh.model.addPhysicalGroup(1, [l_nmc_left, l_lps_left], PG["block_left"],      "block_left")

    gmsh.model.addPhysicalGroup(2, [s_nmc], PG["block_NMC"], "block_NMC")
    gmsh.model.addPhysicalGroup(2, [s_lps], PG["block_LPS"], "block_LPS")

    # -----------------------------------------------------------------
    # Mesh options.  Unstructured quad-dominant mesh - transfinite
    # meshing is awkward on the L-shaped LPS, and the quad-dominant
    # algorithm produces a good-quality mesh that still recombines into
    # nearly all quads.  CRITICALLY: because the arc is one geometric
    # curve shared by both surfaces, the meshes are conformal across
    # it (BreakMeshByBlockGenerator requires this).
    # -----------------------------------------------------------------
    gmsh.option.setNumber("Mesh.Algorithm",            8)   # Frontal-Delaunay for Quads
    gmsh.option.setNumber("Mesh.RecombineAll",         1)   # combine triangles into quads
    gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 3) # Blossom full-quad
    gmsh.option.setNumber("Mesh.ElementOrder",         1)
    gmsh.option.setNumber("Mesh.MshFileVersion",       4.1)

    gmsh.model.mesh.generate(2)
    gmsh.write(filename)

    # -----------------------------------------------------------------
    # Mesh statistics
    # -----------------------------------------------------------------
    node_tags, _, _ = gmsh.model.mesh.getNodes()
    n_quads = sum(
        len(tags)
        for etype, _, tags in zip(
            *(gmsh.model.mesh.getElements(2)[i] for i in range(3))
        )
        if etype == 3   # 4-node quadrangle
    ) if False else None  # fall through to the simpler call below

    # Simpler stats
    elem_types, elem_tags, _ = gmsh.model.mesh.getElements(2)
    n_elems = sum(len(t) for t in elem_tags)

    print(f"[make_mesh] wrote '{filename}'")
    print(f"[make_mesh]   geometry  : {L*1e3:.2f} mm square; NMC arc R = {R*1e3:.4f} mm")
    print(f"[make_mesh]   target h  : {H_EL:.2e} m")
    print(f"[make_mesh]   nodes     : {len(node_tags)}")
    print(f"[make_mesh]   2D elems  : {n_elems}")

    if gui:
        gmsh.fltk.run()

    gmsh.finalize()


if __name__ == "__main__":
    build_mesh(gui=("--gui" in sys.argv or "-gui" in sys.argv))