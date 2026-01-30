# Development Memo for TECM_TEST Enhancement

## Critical Guidelines - ALWAYS FOLLOW

### File Creation Policy
**NEVER create new files without explicit approval from the user.**

**Process for new files:**
1. Identify the need for a new file
2. Ask the user: "I need to create a new file [filename/purpose]. How should I structure/organize this new file?"
3. Wait for user's guidance on file structure, location, and organization
4. Only proceed after receiving explicit instructions

### Collaboration Approach
- Work step by step, case by case
- Preserve existing file structure - do not disturb current organization
- Always ask for guidance on file placement and structure before creating anything new
- This applies to ALL file types: source code, headers, examples, documentation, configuration files, etc.

### Mathematical Equation Display
**ALWAYS display mathematical equations in readable format, NOT LaTeX format.**

**Examples:**
- ✅ GOOD: J_i = -D_i ∇c_i - (z_i F D_i)/(RT) c_i ∇Φ
- ❌ AVOID: $$\mathbf{J}_i = -D_i \nabla c_i - \frac{z_i F D_i}{RT} c_i \nabla \Phi$$

**Guidelines:**
- Use standard mathematical notation with symbols like ∇, ∞, ≈, ≤, ≥
- Use subscripts and superscripts in plain text format (e.g., c_i, x^2)
- Use parentheses for fractions: (numerator)/(denominator)
- Keep equations readable and easy to follow in plain text

### File Testing and Debugging Policy
**ALWAYS test generated files immediately after creation and fix any issues.**

**For MOOSE input files (.i):**
- Run: `./tecm_test-opt -i filename.i --check-input` (syntax check)
- Run: `./tecm_test-opt -i filename.i` (full execution test)
- Debug and fix any parsing errors, missing kernels, or convergence issues
- Ensure the simulation completes successfully

**For LaTeX files (.tex):**
- Run: `pdflatex filename.tex` (compile to PDF)
- Fix any compilation errors or warnings
- Ensure clean PDF generation

**General principle:** Never deliver untested files to the user.

### Reminder Protocol
- Reference this memo in every session
- Remind the user if they forget this guideline
- Always pause and ask before file creation, even if the need seems obvious
- Always use readable mathematical notation, not LaTeX format

## Current Status
- User wants to enhance TECM_TEST framework capabilities
- Focus on making it more useful while respecting existing architecture
- Incremental improvements with user guidance

---
**This memo serves as a persistent reminder of the agreed development workflow.**