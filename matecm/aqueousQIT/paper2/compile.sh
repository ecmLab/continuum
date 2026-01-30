#!/bin/bash
# Compilation script for Paper 2
# Follows AGENTS.md rules: all build artifacts go to build/ folder

set -e  # Exit on error

echo "=========================================="
echo "Compiling Paper 2: Microsphere Experiments"
echo "=========================================="

# Create build directory if it doesn't exist
mkdir -p build

# Clean previous build artifacts
echo "Cleaning previous build..."
rm -f build/*

# Compile with pdflatex (first pass)
echo "Running pdflatex (pass 1)..."
pdflatex -interaction=nonstopmode -output-directory=build paper2_microspheres.tex

# Run bibtex for bibliography
echo "Running bibtex..."
cd build
bibtex paper2_microspheres
cd ..

# Compile with pdflatex (second pass for references)
echo "Running pdflatex (pass 2)..."
pdflatex -interaction=nonstopmode -output-directory=build paper2_microspheres.tex

# Compile with pdflatex (third pass for cross-references)
echo "Running pdflatex (pass 3)..."
pdflatex -interaction=nonstopmode -output-directory=build paper2_microspheres.tex

# Move final PDF to main directory
echo "Moving PDF to main directory..."
mv build/paper2_microspheres.pdf ./

echo ""
echo "=========================================="
echo "Compilation complete!"
echo "Output: paper2_microspheres.pdf"
echo "Build artifacts: build/"
echo "=========================================="
echo ""
echo "To clean build artifacts: rm -rf build/*"
