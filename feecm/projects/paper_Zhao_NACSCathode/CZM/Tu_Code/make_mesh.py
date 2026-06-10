#!/usr/bin/env python3
r"""
make_mesh.py
============

Generate the 2-D NVP / LPS half-cell mesh used by both MOOSE inputs:

  * contact_CZM_3.i  -> cohesive-zone model across the arc.
                        Needs a CONFORMAL interface (shared nodes) so
                        BreakMeshByBlockGenerator can split it.
                        ->  use the default (conformal) mode.

  * contact.i        -> genuine node-to-segment mechanical contact.
                        Needs a NON-CONFORMAL interface: two separate
                        bodies whose arcs are coincident but have
                        independent nodes, so they can actually separate.
                        ->  use the --contact mode.

Geometry:

      (0,L) +------------------------------+ (L,L)
            |                              |
            |                              |
            |          block_NACS           |
            |          (L-shape)           |
            |                              |
       (0,R)*\                             |
            | \.  arc_mid  (the interface, radius R)
            |   \.                         |
            |     \.                       |
            |  block \.                    |
            |  _NVP    \.                  |
            |           \.                 |
            +-------------*----------------+
        (0,0)            (R,0)           (L,0)

    L = 5.0 um                       (full domain side length)
    R = 3.8818230 um                 (NVP quarter-circle radius)

Why the mesh looks the way it does
----------------------------------
ALL of the physics both inputs report lives on the arc (contact pressure,
or cohesive traction / jump / damage, plus the NodalValueSampler that walks
the interface). The NVP disc swells ~2.6% via the thermal/chemical eigenstrain
and is pressed into the soft, yielding LPS under 5 MPa stack pressure, so the
interface is where gradients are sharp and fidelity matters most.

So instead of a uniform mesh (was ~24k elems) we use:

  1. A GRADED SIZE FIELD anchored on the arc. Elements are H_IFACE within a
     REFINE_DIST halo of the interface and coarsen to H_FAR in the bulk. The
     arc lies at zero distance from itself, so it is meshed uniformly at
     H_IFACE -- an even interface discretization for the CZM `sort_by = x`
     sampler and for matched contact surfaces -- without needing a transfinite
     constraint (which is incompatible with full-quad recombination).

  2. H_IFACE is tied to the CZM `interface_thickness` (0.02 um) so one
     element edge ~ the effective cohesive layer thickness.

Conformal vs contact (the only topological difference is the arc)
-----------------------------------------------------------------
  CONFORMAL (default, for CZM): NVP and LPS SHARE the interface arc (one
  curve, opposite orientation), so they share nodes. removeDuplicateNodes()
  welds everything. BreakMeshByBlockGenerator then duplicates the interface
  nodes itself and inserts the cohesive elements.

  CONTACT (--contact, for contact.i): the two surfaces use SEPARATE arcs at
  the same radius, with their own endpoints on the x/y axes. The result is
  two independent bodies that merely touch along the arc and both rest on the
  fixed bottom -- free to separate and slide. removeDuplicateNodes() is NOT
  called, so the coincident interface (and the (R,0)/(0,R) junction) nodes
  stay distinct. block_NVP_right and block_NACS_left land on the two
  independent arcs -> valid primary/secondary contact surfaces.

Sidesets / blocks produced (names match both inputs):
    block_NVP_right, block_NACS_left            (CONTACT mode only; two arcs)
    block_bottom, block_right, block_top, block_left
    block_NVP, block_NACS                       (element blocks)

In CONFORMAL mode the shared arc is left untagged: libMesh's Gmsh reader
forbids one entity from carrying two boundary IDs, and the CZM input gets its
interface sidesets from BreakMeshByBlockGenerator anyway.

Usage
-----
    python make_mesh.py                     # conformal -> CZM filename
    python make_mesh.py --contact           # non-conformal -> contact filename
    python make_mesh.py --out my.msh        # override output name
    python make_mesh.py --gui               # open the gmsh viewer
"""

import sys
import math
import gmsh


# ---------------------------------------------------------------------------
# Geometry / discretization parameters  (all lengths in um (unit system: um-MPa-uN))
# ---------------------------------------------------------------------------
L           = 5.0          # square domain side length
R           = 3.98942280   # NVP quarter-circle radius
H_IFACE     = 0.02         # element edge ON the arc (== CZM interface_thickness)
H_FAR       = 0.20         # target element edge in the bulk far field
REFINE_DIST = 0.10         # halo width around the arc kept fully refined [um]

OUT_FILE_CZM     = "mesh_czm_5050.msh"   # contact_CZM_3.i
OUT_FILE_CONTACT = "mesh_5050.msh"       # contact.i

PG = {
    "block_NVP_right": 1,
    "block_NACS_left":  2,
    "block_bottom":    3,
    "block_right":     4,
    "block_top":       5,
    "block_left":      6,
    "block_NVP":       7,
    "block_NACS":       8,
}


def build_mesh(filename: str | None = None,
               gui: bool = False,
               conformal: bool = True,
               L: float = L,
               R: float = R,
               h_iface: float = H_IFACE,
               h_far: float = H_FAR,
               refine_dist: float = REFINE_DIST) -> None:

    if filename is None:
        filename = OUT_FILE_CZM if conformal else OUT_FILE_CONTACT

    arc_len = 0.5 * math.pi * R                        # quarter-circle length
    n_arc   = max(8, round(arc_len / h_iface)) + 1     # transfinite nodes on arc

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("nmc_lps_contact" if not conformal else "nmc_lps_czm")

    geo = gmsh.model.geo

    # -----------------------------------------------------------------
    # Points common to both modes (origin is the arc centre)
    # -----------------------------------------------------------------
    p_o  = geo.addPoint(0.0, 0.0, 0.0, h_iface)        # origin / arc centre
    p_BR = geo.addPoint(L,   0.0, 0.0, h_far)          # bottom-right corner
    p_TR = geo.addPoint(L,   L,   0.0, h_far)          # top-right corner
    p_TL = geo.addPoint(0.0, L,   0.0, h_far)          # top-left corner

    # -----------------------------------------------------------------
    # Interface points / arcs.
    #   conformal : one shared arc_mid (shared nodes across the interface)
    #   contact   : two coincident arcs with INDEPENDENT endpoints/nodes
    # -----------------------------------------------------------------
    if conformal:
        p_xm_n = p_xm_l = geo.addPoint(R,   0.0, 0.0)  # (R,0) shared
        p_ym_n = p_ym_l = geo.addPoint(0.0, R,   0.0)  # (0,R) shared
        arc_mid_nmc = geo.addCircleArc(p_xm_n, p_o, p_ym_n)
        arc_mid_lps = arc_mid_nmc                       # SAME curve -> conformal
    else:
        p_xm_n = geo.addPoint(R,   0.0, 0.0)            # (R,0) NVP body
        p_xm_l = geo.addPoint(R,   0.0, 0.0)            # (R,0) LPS body (distinct)
        p_ym_n = geo.addPoint(0.0, R,   0.0)            # (0,R) NVP body
        p_ym_l = geo.addPoint(0.0, R,   0.0)            # (0,R) LPS body (distinct)
        arc_mid_nmc = geo.addCircleArc(p_xm_n, p_o, p_ym_n)
        arc_mid_lps = geo.addCircleArc(p_xm_l, p_o, p_ym_l)

    # -----------------------------------------------------------------
    # Curves
    # -----------------------------------------------------------------
    l_nmc_bot = geo.addLine(p_o,    p_xm_n)            # NVP bottom  (origin -> R,0)
    l_lps_bot = geo.addLine(p_xm_l, p_BR)              # LPS bottom  (R,0  -> L,0)
    l_right   = geo.addLine(p_BR,   p_TR)              # right side
    l_top     = geo.addLine(p_TR,   p_TL)             # top side
    l_lps_lft = geo.addLine(p_TL,   p_ym_l)            # LPS left    (0,L -> 0,R)
    l_nmc_lft = geo.addLine(p_ym_n, p_o)               # NVP left    (0,R -> origin)

    # -----------------------------------------------------------------
    # Surfaces
    # -----------------------------------------------------------------
    # NVP quarter disc (pie): bottom, arc, left
    cl_nmc = geo.addCurveLoop([l_nmc_bot, arc_mid_nmc, l_nmc_lft])
    s_nmc  = geo.addPlaneSurface([cl_nmc])

    # LPS L-shape: bottom, right, top, left, then the arc traversed in reverse
    cl_lps = geo.addCurveLoop([l_lps_bot, l_right, l_top, l_lps_lft, -arc_mid_lps])
    s_lps  = geo.addPlaneSurface([cl_lps])

    # Recombine both surfaces to quads (full-quad recombination below
    # guarantees an all-QUAD4 mesh, which Exodus requires per subdomain).
    geo.mesh.setRecombine(2, s_nmc)
    geo.mesh.setRecombine(2, s_lps)

    geo.synchronize()

    # -----------------------------------------------------------------
    # Physical groups  (names consumed by the MOOSE inputs)
    # -----------------------------------------------------------------
    # Interface sidesets:
    #   contact   : two independent arcs -> two valid, separate contact
    #               surfaces (primary block_NVP_right / secondary block_NACS_left).
    #   conformal : the arc is a SINGLE shared curve. libMesh's Gmsh reader
    #               forbids one entity from carrying two boundary IDs, so we do
    #               NOT tag it. The CZM input gets its interface sidesets from
    #               BreakMeshByBlockGenerator (block_NVP_block_NACS etc.) anyway.
    if not conformal:
        gmsh.model.addPhysicalGroup(1, [arc_mid_nmc], PG["block_NVP_right"], "block_NVP_right")
        gmsh.model.addPhysicalGroup(1, [arc_mid_lps], PG["block_NACS_left"],  "block_NACS_left")

    gmsh.model.addPhysicalGroup(1, [l_nmc_bot, l_lps_bot], PG["block_bottom"], "block_bottom")
    gmsh.model.addPhysicalGroup(1, [l_right],              PG["block_right"],  "block_right")
    gmsh.model.addPhysicalGroup(1, [l_top],                PG["block_top"],    "block_top")
    gmsh.model.addPhysicalGroup(1, [l_nmc_lft, l_lps_lft], PG["block_left"],   "block_left")

    gmsh.model.addPhysicalGroup(2, [s_nmc], PG["block_NVP"], "block_NVP")
    gmsh.model.addPhysicalGroup(2, [s_lps], PG["block_NACS"], "block_NACS")

    # -----------------------------------------------------------------
    # Size field: refine on the arc, coarsen into the bulk.
    # -----------------------------------------------------------------
    f_dist = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(f_dist, "CurvesList", [arc_mid_nmc])
    gmsh.model.mesh.field.setNumber(f_dist, "Sampling", 400)

    f_thr = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(f_thr, "InField", f_dist)
    gmsh.model.mesh.field.setNumber(f_thr, "SizeMin", h_iface)
    gmsh.model.mesh.field.setNumber(f_thr, "SizeMax", h_far)
    gmsh.model.mesh.field.setNumber(f_thr, "DistMin", refine_dist)   # fine halo
    gmsh.model.mesh.field.setNumber(f_thr, "DistMax", 0.4 * L)       # full coarsening
    gmsh.model.mesh.field.setAsBackgroundMesh(f_thr)

    # Let the field fully drive sizing.
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints",         0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature",      0)

    # Quad meshing
    gmsh.option.setNumber("Mesh.Algorithm",              8)   # Frontal-Delaunay (quads)
    gmsh.option.setNumber("Mesh.RecombineAll",           1)
    gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 3)   # Blossom full-quad (all-QUAD4)
    gmsh.option.setNumber("Mesh.ElementOrder",           1)
    gmsh.option.setNumber("Mesh.MshFileVersion",         4.1)

    gmsh.model.mesh.generate(2)

    # CRITICAL: weld coincident nodes ONLY for the conformal (CZM) mesh.
    # For contact we must keep the two bodies' interface nodes distinct.
    if conformal:
        gmsh.model.mesh.removeDuplicateNodes()

    gmsh.write(filename)

    # -----------------------------------------------------------------
    # Stats
    # -----------------------------------------------------------------
    node_tags, _, _ = gmsh.model.mesh.getNodes()
    _, elem_tags, _ = gmsh.model.mesh.getElements(2)
    n_elems  = sum(len(t) for t in elem_tags)
    frac_nmc = (0.25 * math.pi * R * R) / (L * L)
    mode     = "CONFORMAL (CZM)" if conformal else "NON-CONFORMAL (contact)"

    print(f"[make_mesh] wrote '{filename}'")
    print(f"[make_mesh]   mode          : {mode}")
    print(f"[make_mesh]   domain        : {L:.2f} mm square")
    print(f"[make_mesh]   NVP arc R     : {R:.6f} mm  (area fraction {frac_nmc:5.1%})")
    print(f"[make_mesh]   interface h   : {h_iface:.3e} mm  ({n_arc-1} elems on arc)")
    print(f"[make_mesh]   refine halo   : {refine_dist:.3f} mm")
    print(f"[make_mesh]   far-field h   : {h_far:.3e} mm")
    print(f"[make_mesh]   nodes         : {len(node_tags)}")
    print(f"[make_mesh]   2D elements   : {n_elems}")
    if not conformal:
        print("[make_mesh]   note          : interface nodes are duplicated; "
              "do NOT run BreakMeshByBlockGenerator on this mesh.")

    if gui:
        gmsh.fltk.run()

    gmsh.finalize()


def _parse_args(argv):
    conformal = "--contact" not in argv
    gui = ("--gui" in argv or "-gui" in argv)
    out = None
    if "--out" in argv:
        out = argv[argv.index("--out") + 1]
    return out, gui, conformal


if __name__ == "__main__":
    out_file, show_gui, is_conformal = _parse_args(sys.argv)
    build_mesh(filename=out_file, gui=show_gui, conformal=is_conformal)