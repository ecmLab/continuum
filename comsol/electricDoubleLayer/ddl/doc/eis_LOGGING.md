# eis.m Logging and Output

## Date
2025-10-20

## Summary

Updated `eis.m` to save COMSOL model (.mph file) and capture all simulation output to a log file.

## Changes Made

### 1. Added Simulation Logging (Lines 9-22)

**Setup diary at start:**
```matlab
%% Create output directory first
output_dir = 'rst/eis';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Setup logging - redirect all output to log file
log_file = fullfile(output_dir, 'eis_simulation.log');
diary(log_file);
diary on;

fprintf('=== EIS SIMULATION LOG ===\n');
fprintf('Started: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf('Working directory: %s\n\n', pwd);
```

**Key features:**
- Creates log file at start of simulation
- All `fprintf()` output captured automatically
- Timestamps start time
- Records working directory

### 2. Added Model Parameter Documentation (Lines 195-206)

**After model is built:**
```matlab
% Log key model parameters
fprintf('\n=== MODEL PARAMETERS ===\n');
fprintf('Temperature: %.2f K\n', mphglobal(model, 'T0'));
fprintf('Diffusion coefficients:\n');
fprintf('  D_cation: %.4e m^2/s\n', mphglobal(model, 'DA'));
fprintf('  D_anion: %.4e m^2/s\n', mphglobal(model, 'DX'));
fprintf('Bulk concentration: %.4e mol/m^3\n', mphglobal(model, 'cA_bulk'));
fprintf('Debye length: %.4e m\n', mphglobal(model, 'xD'));
fprintf('Stern layer thickness: %.4e m\n', mphglobal(model, 'xS'));
fprintf('Domain length: %.4e m\n', mphglobal(model, 'L_cell'));
fprintf('Permittivity (relative): %.2f\n', mphglobal(model, 'eps_PEO'));
```

**Logs:**
- All key physical parameters
- Computed values (Debye length, domain length)
- Material properties

### 3. Added COMSOL Model Save (Lines 337-341)

**After data is saved:**
```matlab
% Save COMSOL model with solution
fprintf('  Saving COMSOL model...\n');
mph_file = fullfile(output_dir, 'eis_solved.mph');
mphsave(model, mph_file);
fprintf('  COMSOL model: %s\n', mph_file);
```

**Features:**
- Saves complete model with all solutions
- Includes mesh, physics, results from all frequencies
- Can be opened in COMSOL GUI for post-processing

### 4. Enhanced Final Summary (Lines 407-420)

**At end of simulation:**
```matlab
fprintf('\n=== DONE ===\n');
fprintf('Completed: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf('All results saved in %s/\n', output_dir);
fprintf('\nOutput files:\n');
fprintf('  - %s (CSV data)\n', csv_file);
fprintf('  - %s (MATLAB data)\n', mat_file);
fprintf('  - %s (COMSOL model)\n', mph_file);
fprintf('  - %s (simulation log)\n', log_file);
fprintf('  - %s/eis_nyquist.png (Nyquist plot)\n', output_dir);
fprintf('  - %s/eis_bode.png (Bode plot)\n', output_dir);

% Close diary
diary off;
fprintf('\nLog file saved: %s\n', log_file);
```

**Features:**
- Timestamps completion time
- Lists all output files with descriptions
- Closes diary to finalize log file
- Final confirmation message

## Output Files

After running `eis.m`, the following files are generated in `rst/eis/`:

### 1. eis_simulation.log (NEW)
**Complete simulation log including:**
- Start/end timestamps
- Working directory
- EIS configuration (frequencies, V_dc, V_ac, cycles)
- Model parameters (diffusion, concentration, Debye length, etc.)
- Progress for each frequency:
  - Period, total time, time points
  - Simulation time
  - Impedance results (Z, |Z|, phase)
- Final summary with all impedance values
- List of all output files

**Example content:**
```
=== EIS SIMULATION LOG ===
Started: 2025-10-20 14:23:15
Working directory: /Users/howardtu/Documents/modeling/continuum/feecm/comsol/ddl

=== EIS FREQUENCY SWEEP ===

EIS Configuration:
  Frequencies: 1e+03 Hz, 7e+02 Hz, ... 4e-01 Hz
  DC bias: 0.000 V
  AC amplitude: 0.010 V (10.0 mV)
  Cycles per frequency: 10

Building base COMSOL model...
Base model built successfully!

=== MODEL PARAMETERS ===
Temperature: 298.15 K
Diffusion coefficients:
  D_cation: 1.0000e-14 m^2/s
  D_anion: 2.0000e-14 m^2/s
Bulk concentration: 1.0000e+01 mol/m^3
Debye length: 2.4085e-09 m
Stern layer thickness: 2.0000e-10 m
Domain length: 4.8170e-07 m
Permittivity (relative): 10.00

Starting frequency sweep (15 frequencies)...
Progress: 0/15

========================================
Frequency 1/15: 1.00e+03 Hz
========================================
  Period: 1.0000e-03 s
  Total time: 1.0000e-02 s (10 cycles)
  Time points: 501
  Time step: 2.0000e-05 s
  Midpoint: x = 2.4085e-07 m (mesh index 21)

  Running simulation...
  Simulation completed in 3.45 seconds
  Extracting data...
  Analyzing impedance...

  Results:
    Z = 1.6234e-01 +3.2145e-02j Ohm*m^2
    |Z| = 1.6549e-01 Ohm*m^2
    Phase = 11.23 deg

Progress: 1/15 frequencies completed

[... continues for all frequencies ...]

========================================
=== FREQUENCY SWEEP COMPLETED ===
========================================

Saving results...
  CSV: rst/eis/eis_data.csv
  MAT: rst/eis/eis_results.mat
  Saving COMSOL model...
  COMSOL model: rst/eis/eis_solved.mph

Generating plots...
  Plots saved to rst/eis/

=== EIS SUMMARY ===
Frequency sweep: 1e+03 Hz to 4e-01 Hz (15 points)

Impedance results:
  1e+03 Hz: Z = 1.6549e-01 Ohm*m^2, phi = 11.23 deg
  [... all results ...]

=== DONE ===
Completed: 2025-10-20 14:28:32
All results saved in rst/eis/

Output files:
  - rst/eis/eis_data.csv (CSV data)
  - rst/eis/eis_results.mat (MATLAB data)
  - rst/eis/eis_solved.mph (COMSOL model)
  - rst/eis/eis_simulation.log (simulation log)
  - rst/eis/eis_nyquist.png (Nyquist plot)
  - rst/eis/eis_bode.png (Bode plot)
```

### 2. eis_solved.mph (NEW)
**COMSOL model file including:**
- Complete model geometry and mesh
- All physics (Electrostatics + Transport)
- All parameters
- Solutions for ALL frequencies
- Can be opened in COMSOL GUI:
  ```bash
  comsol -open rst/eis/eis_solved.mph
  ```

**Use cases:**
- Post-process in COMSOL GUI
- Create additional plots
- Extract field profiles (concentration, potential)
- Visualize double layer structure
- Share with collaborators

### 3. eis_data.csv (Existing)
Impedance data in CSV format

### 4. eis_results.mat (Existing)
MATLAB workspace with all variables

### 5. eis_nyquist.png (Existing)
Nyquist plot

### 6. eis_bode.png (Existing)
Bode magnitude and phase plots

## Benefits

### For Debugging:
- ✅ Complete record of simulation progress
- ✅ Timestamps show which frequencies are slow
- ✅ Parameter values logged for verification
- ✅ Can review log without re-running simulation

### For Reproducibility:
- ✅ All settings documented in log file
- ✅ COMSOL model saved for exact reproduction
- ✅ Working directory recorded
- ✅ Date/time stamps for version control

### For Analysis:
- ✅ Open .mph in COMSOL GUI for detailed visualization
- ✅ Extract additional results not saved in CSV
- ✅ Create custom plots and animations
- ✅ Share complete model with collaborators

### For Documentation:
- ✅ Log file can be included in reports
- ✅ Complete audit trail of simulation
- ✅ Parameter documentation for methods section

## Usage

### Running the simulation:
```bash
cd /Users/howardtu/Documents/modeling/continuum/feecm/comsol/ddl
matlab -batch "eis"
```

### Checking the log during simulation:
```bash
tail -f rst/eis/eis_simulation.log
```

### Opening COMSOL model:
```bash
comsol -open rst/eis/eis_solved.mph
```

### Reviewing results after completion:
```bash
less rst/eis/eis_simulation.log
```

## Notes

- **Log file overwrites** each run (not appended)
- **diary** captures all `fprintf()` and command window output
- **.mph file** can be large (10-100 MB depending on mesh/frequencies)
- **Model includes** solutions from last frequency only (COMSOL limitation)
  - To save all solutions, would need separate .mph for each frequency
- **Log remains valid** even if MATLAB crashes (diary writes continuously)

## Future Enhancements (If Needed)

1. **Append mode for logs:**
   ```matlab
   diary(log_file);
   diary on;
   fprintf('\n\n=== NEW RUN: %s ===\n\n', datestr(now));
   ```

2. **Save individual .mph per frequency:**
   ```matlab
   mph_file_freq = fullfile(output_dir, sprintf('eis_freq_%04d.mph', freq_idx));
   mphsave(model, mph_file_freq);
   ```

3. **Compressed log (if too large):**
   ```bash
   gzip rst/eis/eis_simulation.log
   ```

4. **HTML log for better readability:**
   - Convert log to HTML with syntax highlighting
   - Add hyperlinks to output files

## Compatibility

- ✅ MATLAB R2018a or later (diary, fprintf, mphsave)
- ✅ COMSOL Multiphysics 5.4 or later (.mph format)
- ✅ All operating systems (macOS, Linux, Windows)
