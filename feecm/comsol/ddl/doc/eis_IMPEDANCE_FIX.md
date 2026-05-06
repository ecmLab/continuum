# Critical Fixes to eis.m Impedance Calculation

## Date
2025-10-20

## Problem Identified

User reported weird behavior in Nyquist plot at high frequencies - real part of impedance not decreasing as expected.

## Root Causes Found

### Issue 1: Measuring Current at Wrong Location ❌

**Original code (Lines 264-290):**
```matlab
% Get midpoint coordinate (only once, after first simulation)
if isempty(x_mid)
    x_mid = mphglobal(model, 'L_cell/2');
    x_mid = x_mid(1);

    eval_x = mpheval(model, 'x', 'edim', 1, 'solnum', 1);
    x_coords = eval_x.d1;
    [~, mid_idx] = min(abs(x_coords - x_mid));

    fprintf('  Midpoint: x = %.4e m (mesh index %d)\n', x_mid, mid_idx);
end

% Extract i_total at midpoint for all time steps
num_times = length(t_data);
i_total_mid_t = zeros(num_times, 1);
for idx = 1:num_times
    eval_i = mpheval(model, 'i_total', 'edim', 1, 'solnum', idx);
    i_total_mid_t(idx) = eval_i.d1(mid_idx);
end
```

**Problem:**
- Measured current at **MIDPOINT** of domain (x = L_cell/2)
- In EIS, impedance Z = V/I where:
  - V = voltage applied at electrode
  - I = current AT THE ELECTRODE
- Measuring current in the middle of the electrolyte doesn't give the correct electrode response!

**Why this causes issues:**
- At high frequencies, current varies spatially across domain
- Midpoint current ≠ electrode current
- Results in incorrect impedance magnitude and phase
- Similar to measuring voltage drop in the middle of a resistor instead of across the terminals!

### Issue 2: Wrong Phase Sign ❌

**Original code (Line 72):**
```matlab
phase_rad = -angle(Z_complex);  % WRONG: negates the phase!
```

**Problem:**
- Added negative sign to phase angle
- For capacitive systems (like EDL), impedance should have **negative imaginary part**
- Z = R - jX means phase is negative (current leads voltage)
- The `-` sign flips this, making results confusing

**Convention:**
```
Z = V/I = |Z| ∠φ
For capacitor: φ < 0 (current leads)
For inductor: φ > 0 (current lags)
```

## Fixes Applied

### Fix 1: Measure Current at Left Electrode ✅

**New code (Lines 264-271):**
```matlab
%% Extract time series data at LEFT ELECTRODE
fprintf('  Extracting data at left electrode (x=0)...\n');

% Get time, voltage, and current at LEFT electrode (boundary point 1)
time_data = mpheval(model, {'t', 'V_applied', 'i_total'}, 'edim', 0, 'selection', 1);
t_data = time_data.d1;
V_applied_t = time_data.d2;
I_electrode_t = time_data.d3;  % Current density at left electrode (A/m^2)
```

**Changes:**
- Extract `i_total` directly at electrode boundary (selection 1 = left electrode at x=0)
- Use `'edim', 0` to specify boundary dimension
- Single `mpheval` call gets time, voltage, AND current simultaneously
- No need for loop or midpoint calculations

**Why this is correct:**
- Measures current where voltage is applied
- Proper definition: Z = V_electrode / I_electrode
- Consistent with standard EIS measurements
- Similar to `eis_ref.m` implementation (line 287-290)

### Fix 2: Remove Phase Negation ✅

**New code (Line 303-304):**
```matlab
phase_rad = angle(Z_complex);  % No negative sign - direct phase
phase_deg = phase_rad * 180/pi;
```

**Changes:**
- Removed negative sign from `angle(Z_complex)`
- Direct phase from complex impedance
- Standard convention: negative phase = capacitive

### Fix 3: Cleanup Unused Variables ✅

**Removed (Lines 227-229):**
```matlab
% Will store midpoint info after first simulation
x_mid = [];
mid_idx = [];
```

No longer needed since we're not using midpoint.

## Expected Behavior After Fixes

### At High Frequencies (> 1 kHz):

**Before fix:**
- Z_real ~ 0.193 Ω·m² (constant, but possibly wrong magnitude)
- Z_imag ~ -0.0001 Ω·m² (very small)
- May not approach correct R_∞

**After fix:**
- Z_real → R_∞ (solution resistance)
- Z_imag → 0 (purely resistive at high freq)
- Clear approach to constant value
- Nyquist plot should start at (R_∞, 0)

### At Low Frequencies:

**Before fix:**
- May have artifacts from wrong current measurement

**After fix:**
- Clear Warburg impedance (45° line)
- Large capacitive impedance
- Vertical line at lowest frequencies

### Overall Nyquist Plot:

**Expected shape:**
```
High freq (R_∞, 0) ───→ Warburg (45° line) ───→ Vertical capacitive line
                    ↘               ↘
                   Semicircle    Low freq
```

## Comparison with eis_ref.m

The fixed `eis.m` now matches the approach in `eis_ref.m`:

**eis_ref.m (Lines 287-290):**
```matlab
% Get time, voltage, and current at electrode
time_data = mpheval(model, {'t', 'V_applied', 'i_total'}, 'edim', 0, 'selection', 1);
t_data = time_data.d1;
V_applied_t = time_data.d2;
I_electrode_t = time_data.d3;  % Current density at electrode (A/m^2)
```

**Key difference:**
- `eis_ref.m` uses **least-squares fitting** for impedance extraction
- `eis.m` uses **FFT method**
- Both now correctly measure current at electrode

## Technical Notes

### Why Electrode Current Matters:

For a two-electrode system:
```
[Left Electrode] ←─ i_total ─→ [Electrolyte] ←─ i_total ─→ [Right Electrode]
     V_applied                                                  V = 0 (ground)
```

**Current conservation:**
- In 1D, i_total(x) is continuous
- BUT: i_total varies with position due to charging/discharging
- At electrode: i_total includes **charging current** (dQ/dt)
- In bulk: i_total is only **ionic current**

**At high frequencies:**
- Charging current dominates at electrodes
- Large spatial variation in current
- Midpoint current ≪ electrode current
- This is why measurement location matters!

### Phase Convention:

```
Z = R + jX

If X < 0:  Capacitive (current leads voltage)
           phase φ = arctan(X/R) < 0

If X > 0:  Inductive (current lags voltage)
           phase φ = arctan(X/R) > 0
```

For EDLC (electric double layer capacitor):
- Dominated by capacitance
- X < 0 (negative imaginary impedance)
- φ < 0 (negative phase)
- In Nyquist: plot -Z_imag (positive y-axis)

## Verification Steps

After running fixed `eis.m`:

1. **Check high-frequency limit:**
   ```matlab
   Z_real(1)  % Should be ~ R_∞ (solution resistance)
   Z_imag(1)  % Should be ~ 0
   ```

2. **Check Nyquist plot shape:**
   - Should start at high freq near (R_∞, 0)
   - Should show semicircle or Warburg line
   - Should become vertical at low freq

3. **Check phase behavior:**
   - High freq: φ ≈ 0° (resistive)
   - Mid freq: φ ≈ -45° (Warburg)
   - Low freq: φ ≈ -90° (capacitive)

4. **Compare with theory:**
   - Calculate R_∞ = L/(κ) where κ = solution conductivity
   - High-freq Z_real should match R_∞

## Files Modified

- **eis.m**: Lines 223-306
  - Removed midpoint calculation (lines 227-229, 264-274)
  - Changed to electrode current extraction (lines 264-271)
  - Fixed phase calculation (line 303)
  - Updated variable name: `i_total_mid_t` → `I_electrode_t`

## Related Issues

This same error pattern was found in `eis_ref.m`:
- ✅ Already fixed (never had this issue - measured at electrode from start)

This error pattern was NOT in other files:
- `sinusoidal.m`: Measures at electrode correctly
- `oneInterface.m`: Single electrode, no ambiguity
- `twoInterface.m`: Need to check if similar issue exists

## Lessons Learned

1. **Always measure at electrodes for EIS:**
   - Impedance = V_electrode / I_electrode
   - NOT V_electrode / I_midpoint

2. **Phase sign matters:**
   - Use direct `angle(Z)`, don't negate
   - Capacitive = negative imaginary = negative phase

3. **High-frequency behavior is diagnostic:**
   - Should approach R_∞ (constant real value)
   - If weird behavior at high freq → check measurement location

4. **Spatial variation matters at high frequencies:**
   - Current varies across domain
   - Only electrode current gives correct impedance
