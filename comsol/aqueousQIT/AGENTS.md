
1) LaTeX Build Hygiene
- When compiling `.tex` documents, only keep the `.tex` source and the resulting `.pdf` in the document's folder.
- Move all other LaTeX build artifacts (e.g., `.aux`, `.log`, `.out`, `.toc`, `.synctex.gz`, `.fls`, `.fdb_latexmk`, `.bbl`, `.blg`) into a subfolder named `build` inside that same folder. Create the folder if it does not exist.
- If multiple `.tex` files share a folder, they may share the same `build/` subfolder for their auxiliary files.
