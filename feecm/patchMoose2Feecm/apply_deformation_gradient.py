#!/usr/bin/env python3
"""
apply_deformation_gradient.py  --  feecm MOOSE patch (robust & idempotent)
===========================================================================
Adds the `<base>deformation_gradient` AD material property to MOOSE's
incremental finite-strain calculators. feecm's swelling / visco materials
(ecm_test/tecm_test) consume it via:

    getADMaterialProperty<RankTwoTensor>(_base_name + "deformation_gradient")

WHY A SCRIPT INSTEAD OF COPYING WHOLE FILES
-------------------------------------------
The old `patchMoose2Feecm/*.C,*.h` approach copied entire OLD-MOOSE source
files over the current ones. That clobbers every unrelated upstream change and
broke the build when MOOSE:
  * renamed the module  tensor_mechanics -> solid_mechanics
  * switched AD math to ADL  (`using std::cbrt; cbrt(...)`  not  `std::cbrt`)
  * added a Jacobian-singularity perturbation for F == I

This script instead makes ONLY the minimal additive insertions, anchored on
stable code landmarks (whitespace-tolerant), and is IDEMPOTENT: re-run it any
time, in particular after every `git pull` / rebuild of MOOSE. If an anchor
can no longer be found, it stops with a clear message rather than silently
producing a broken tree.

USAGE
-----
    python3 apply_deformation_gradient.py            # auto-detect MOOSE, apply
    python3 apply_deformation_gradient.py --check     # report only, no writes
    python3 apply_deformation_gradient.py --moose-dir /path/to/moose
    # revert (restore pristine MOOSE) is just:  git -C $MOOSE_DIR checkout -- modules/solid_mechanics
"""
import argparse, os, re, sys

# module-relative paths of the files we touch (stable across recent MOOSE)
FILES = [
    ("base_h", "include/materials/ADComputeIncrementalStrainBase.h"),
    ("base_c", "src/materials/ADComputeIncrementalStrainBase.C"),
    ("strain", "src/materials/ADCompute1DFiniteStrain.C"),
    ("strain", "src/materials/ADCompute2DFiniteStrain.C"),
    ("strain", "src/materials/ADComputeFiniteStrain.C"),
]


def detect_moose_dir(explicit):
    if explicit:
        return explicit
    if os.environ.get("MOOSE_DIR"):
        return os.environ["MOOSE_DIR"]
    here = os.path.dirname(os.path.abspath(__file__))
    mk = os.path.join(here, "..", "ecm_test", "Makefile")
    if os.path.isfile(mk):
        with open(mk) as f:
            for line in f:
                m = re.match(r"\s*MOOSE_DIR\s*:?=\s*(\S+)", line)
                if m:
                    return m.group(1)
    sys.exit("ERROR: cannot locate MOOSE. Set $MOOSE_DIR or pass --moose-dir.")


def find_module_dir(moose_dir):
    for name in ("solid_mechanics", "tensor_mechanics"):  # new name, then legacy
        d = os.path.join(moose_dir, "modules", name)
        if os.path.isdir(d):
            return d
    sys.exit(f"ERROR: no solid_mechanics/tensor_mechanics module under {moose_dir}/modules")


# ---- transforms: each takes the file text, returns (new_text, message) ----
#      They are individually idempotent (guarded by presence checks).

def patch_base_h(t):
    msgs = []
    # (1) add `_deformation_gradient` to the `using` member-import macro so the
    #     templated 3D calculator can see it (two-phase name lookup).
    if "::_deformation_gradient;" not in t:
        m = re.search(r"^( *using +(\w+<R2>)::_rotation_increment;)( *)\\$", t, re.M)
        if not m:
            raise RuntimeError("anchor not found: macro `using ...::_rotation_increment; \\`")
        full, cls = m.group(0), m.group(2)
        line = f"  using {cls}::_deformation_gradient;"
        line = line + " " * max(1, len(full) - 1 - len(line)) + "\\"  # align trailing backslash
        t = t[:m.end()] + "\n" + line + t[m.end():]
        msgs.append("macro using-line")
    # (2) declare the protected member next to _rotation_increment
    if re.search(r"&\s*_deformation_gradient\s*;", t) is None:
        m = re.search(r"^( *)ADMaterialProperty<RankTwoTensor> & _rotation_increment;$", t, re.M)
        if not m:
            raise RuntimeError("anchor not found: member `ADMaterialProperty<RankTwoTensor> & _rotation_increment;`")
        ind = m.group(1)
        t = t[:m.end()] + f"\n\n{ind}ADMaterialProperty<RankTwoTensor> & _deformation_gradient;" + t[m.end():]
        msgs.append("member decl")
    return t, (", ".join(msgs) if msgs else "already present")


def patch_base_c(t):
    msgs = []
    # (1) declare the AD property in the constructor initializer list
    if '"deformation_gradient"' not in t:
        anchor = '"rotation_increment")),\n'
        i = t.find(anchor)
        if i < 0:
            raise RuntimeError('anchor not found: `"rotation_increment")),`')
        j = i + len(anchor)
        t = (t[:j]
             + '    _deformation_gradient(\n'
               '        this->template declareADProperty<RankTwoTensor>(_base_name + "deformation_gradient")),\n'
             + t[j:])
        msgs.append("ctor declare")
    # (2) initialise it to identity in initQpStatefulProperties()
    if "_deformation_gradient[_qp].setToIdentity()" not in t:
        anchor = "  _total_strain[_qp].zero();\n"
        i = t.find(anchor)
        if i < 0:
            raise RuntimeError("anchor not found: `_total_strain[_qp].zero();` in initQpStatefulProperties")
        j = i + len(anchor)
        t = t[:j] + "  _deformation_gradient[_qp].setToIdentity();\n" + t[j:]
        msgs.append("initQp identity")
    return t, (", ".join(msgs) if msgs else "already present")


def patch_strain(t):
    # Store F = I + grad(u) at the Gauss point, just before `A -= Fbar` (where A
    # still holds the *current* displacement gradient). Works for 1D/2D/3D.
    if "_deformation_gradient[_qp] = A;" in t:
        return t, "already present"
    m = re.search(r"^([ \t]*)A -= Fbar;", t, re.M)
    if not m:
        raise RuntimeError("anchor not found: `A -= Fbar;`")
    ind = m.group(1)
    block = (f"{ind}_deformation_gradient[_qp] = A;\n"
             f"{ind}_deformation_gradient[_qp].addIa(1.0); // Gauss point deformation gradient\n\n")
    t = t[:m.start()] + block + t[m.start():]
    return t, "inserted F = I + grad(u)"


HANDLERS = {"base_h": patch_base_h, "base_c": patch_base_c, "strain": patch_strain}


def main():
    ap = argparse.ArgumentParser(description="Apply feecm _deformation_gradient patch to MOOSE.")
    ap.add_argument("--moose-dir", help="MOOSE root (default: $MOOSE_DIR or ecm_test/Makefile)")
    ap.add_argument("--check", action="store_true", help="report status only; do not modify files")
    args = ap.parse_args()

    moose = detect_moose_dir(args.moose_dir)
    mod = find_module_dir(moose)
    print(f"MOOSE   : {moose}")
    print(f"module  : {mod}")
    print(f"mode    : {'CHECK (no writes)' if args.check else 'APPLY'}\n")

    rc = 0
    for kind, rel in FILES:
        path = os.path.join(mod, rel)
        name = os.path.basename(rel)
        if not os.path.isfile(path):
            print(f"  MISSING      {name}")
            rc = 1
            continue
        with open(path) as f:
            src = f.read()
        try:
            new, msg = HANDLERS[kind](src)
        except RuntimeError as e:
            print(f"  ANCHOR-FAIL  {name}: {e}")
            rc = 1
            continue
        changed = new != src
        if args.check:
            print(f"  {'NEEDS-PATCH' if changed else 'ok         '}  {name}  ({msg})")
        else:
            if changed:
                with open(path, "w") as f:
                    f.write(new)
            print(f"  {'PATCHED    ' if changed else 'ok         '}  {name}  ({msg})")

    if rc:
        print("\nFAILED: one or more files could not be patched.")
        print("MOOSE's structure may have changed; review the anchors above and update this script.")
    else:
        print("\nCHECK complete." if args.check else "\nAll patches applied. Rebuild with `make -j` in the app dir.")
    sys.exit(rc)


if __name__ == "__main__":
    main()
