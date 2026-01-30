# Agent Guidelines (Repo-Wide)

These rules apply to all coding agents working in this repository.
Scope: entire repo.

## Core Rules for plan and strategies

1) Follow the Roadmap
- Use the project plan at `doc/codex/TECM_Project_Plan.md` to guide priorities and phase scope. Link your work to the relevant phase and acceptance criteria from the project plan
- Keep a repo‑wide technical document `doc/codex/TECM_Development.tex` (PDF built alongside) for the overall framework.
- Communication: When uncertain about scope, naming, or placement, ask before proceeding.

2) Phase‑Oriented Development Layout
- Implement work phase by phase. Avoid opportunistic refactors outside the active phase.
- For each phase N, add:
  - A phase technical file `doc/codex/phaseN_<name>.tex` (e.g., `phase2_edl.tex`) and build its `phaseN_<name>.pdf` in the same folder (apply Rule 6 for build artifacts).
  - Phase‑scoped code folders: `src/phaseN_<name>/`, `include/phaseN_<name>/`, and inputs/examples in `examples/phaseN_<name>/`.
  - Within each phase code folder, place files into appropriate subfolders by type: `kernels/` for Kernels, `bcs/` for boundary conditions, `materials/` for materials, `interfacekernels/` for interface kernels, etc. Example: `src/phase2_edl/bcs/SternRobinBC.C`, `include/phase2_edl/bcs/SternRobinBC.h`.
- Use consistent lowercase names (e.g., `phase2_edl`), and mirror this naming across doc/src/include/examples.
- When starting a new phase, copy the prior phase structure and update: governing equations, MOOSE weak forms, coding plan, and references to acceptance criteria in the project plan.

## Notes & Conventions for coding environment
3) Verify What You Create (No untested hand‑offs)
- Any input or source file you add or modify must be checked with the appropriate executable command(s) before presenting it to the user.
- Minimum checks (examples):
  - MOOSE input files (`.i`):
    - `../../tecm_test-opt -i <file>.i --check-input`
    - `../../tecm_test-opt -i <file>.i`
  - LaTeX (`.tex`):
    - `pdflatex <file>.tex` (twice if references) and confirm the PDF is generated
  - C++ sources/headers:
    - Build with the repo’s build system (e.g., `make -j`), fix errors/warnings relevant to your change
  - Scripts:
    - Ensure executable bit and a minimal run (`--help` or a smoke test)
- If a file cannot be executed/compiled in this environment, state why and provide a concrete, reproducible path for the user to validate locally.

4) Use the MOOSE Conda Environment
- Before compiling/building the MOOSE executable (e.g., `make`, `make -j4`), ensure the `moose` conda environment is active.
- Example activation: `source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate moose` or run commands via `conda run -n moose ...`.

5) Run Executables From the Input Directory
- When running the executable with an input file, first `cd` to the directory containing the `.i` file.
- Invoke the executable by its correct location relative to that directory. For example: `../../tecm_test-opt -i <file>.i --check-input` followed by `../../tecm_test-opt -i <file>.i`.
- Adjust the relative path as needed, but always run from the input file’s directory and point to the correct executable path.

6) LaTeX Build Hygiene
- When compiling `.tex` documents, only keep the `.tex` source and the resulting `.pdf` in the document’s folder.
- Move all other LaTeX build artifacts (e.g., `.aux`, `.log`, `.out`, `.toc`, `.synctex.gz`, `.fls`, `.fdb_latexmk`, `.bbl`, `.blg`) into a subfolder named `build` inside that same folder. Create the folder if it does not exist.
- If multiple `.tex` files share a folder, they may share the same `build/` subfolder for their auxiliary files.

7) Implementing MOOSE Boundary Conditions
- Do not implement BCs like element kernels. Kernels assemble volumetric weak forms; BC objects supply boundary residuals directly on the boundary.
- For diffusion/Poisson, the natural boundary term from integration by parts is handled exclusively by BC objects. Replace the eliminated flux using the BC relation and add that residual on the boundary. Do not add boundary integrals inside volumetric kernels.
- Example (Stern/Robin at an electrode): with normal pointing out of the electrolyte, if `epsilon * ∂n φ = C_S (φ_M − φ)`, the BC should add `− ∫ C_S (φ_M − φ) w dA` on that boundary.
- Prefer AD variants (e.g., derive from `ADIntegratedBC`) and pass physical parameters (e.g., `epsilon0`, `stern_relative_permittivity`, `stern_thickness`, `metal_potential`); compute derived quantities (e.g., `C_S = epsilon0 * stern_relative_permittivity / stern_thickness`) inside the BC.
