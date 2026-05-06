
1) LaTeX Build Hygiene
- When compiling `.tex` documents, only keep the `.tex` source and the resulting `.pdf` in the document's folder.
- Move all other LaTeX build artifacts (e.g., `.aux`, `.log`, `.out`, `.toc`, `.synctex.gz`, `.fls`, `.fdb_latexmk`, `.bbl`, `.blg`) into a subfolder named `build` inside that same folder. Create the folder if it does not exist.
- If multiple `.tex` files share a folder, they may share the same `build/` subfolder for their auxiliary files.

2) Folder and File Organization
- **Main MATLAB files**: Keep one `.m` file per development phase in the root directory (`ddl/`)
  - Phase 0: `ddl0.m` (reference, do not modify)
  - Phase 1: `oneInterface.m` (renamed from `ddl_oneInterface.m`)
  - Phase 2: `twoInterface.m` (renamed from `ddl_twoInterface.m`)
  - Phase 3: `sinusoidal.m` (sinusoidal AC drive study)
  - Phase 4: `eis.m` (EIS analysis)
  - Phase 5: `pulse.m` (future development)

- **Documentation**: All documentation files (`.md`, `.tex`, `.pdf`) go in `doc/` folder

- **Results**: All simulation results organized by phase in `rst/` with subfolders:
  - `rst/oneInterface/` - Phase 1 results
  - `rst/twoInterface/` - Phase 2 results
  - `rst/sinusoidal/` - Phase 3 results
  - `rst/eis/` - Phase 4 results
  - `rst/pulse/` - Phase 5 results (when ready)

- **File Naming Convention**: All output files follow the pattern `{phaseName}_{description}.{ext}`
  - Model files: `{phase}_solved.mph`
  - Results: `{phase}_results.mat`
  - Data: `{phase}_data*.csv`
  - Plots: `{phase}_{plottype}.png`

- **Development Workflow**:
  - Only one `.m` file per phase - keep editing that single file until phase is complete
  - Each phase writes outputs to its own `rst/{phaseName}/` directory
  - Do not create multiple versions of the same phase file
