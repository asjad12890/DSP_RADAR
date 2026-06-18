# 📡 Radar DSP Processing Framework
**CE363 — Digital Signal Processing | Complex Engineering Problem**
*M. Asjad | 2023386 | GIKI*

---

## Overview

A complete radar signal processing simulation built in MATLAB, covering the full pipeline from waveform generation to target detection. The system models three radar waveforms, four CFAR detection schemes, three clutter mitigation filters, and five Swerling target fluctuation models — with parametric analyses comparing their performance across SNR, bandwidth, and CFAR parameters.

---

## Project Structure

```
DSP_RADAR/
│
├── main_radar.m              ← Entry point — interactive simulation runner
│
├── Waveform Generation
│   ├── generate_LFM.m        ← Linear Frequency Modulated (LFM) chirp
│   ├── generate_PhaseCode.m  ← Barker-13 phase-coded waveform
│   └── generate_Hybrid.m     ← LFM modulated by Barker-13 phase code
│
├── Clutter Mitigation
│   ├── apply_MTI.m           ← 3-pulse canceller (Moving Target Indicator)
│   ├── apply_adaptive.m      ← LMS adaptive filter (clutter subtraction)
│   └── apply_doppler.m       ← Doppler notch filter (zero-Doppler rejection)
│
├── Detection
│   ├── detect_fixed.m        ← Fixed threshold detector (baseline)
│   ├── detect_CA_CFAR.m      ← Cell-Averaging CFAR
│   ├── detect_OS_CFAR.m      ← Ordered-Statistics CFAR
│   └── detect_GO_CFAR.m      ← Greatest-Of CFAR
│
└── Analysis Scripts
    ├── plot_ambiguity.m          ← 2D ambiguity functions for all waveforms
    ├── analysis_SNR_vs_Pd.m      ← SNR vs Probability of Detection (Monte Carlo)
    ├── analysis_bandwidth.m      ← Bandwidth vs range resolution and sidelobes
    ├── analysis_CFAR_params.m    ← CFAR training cells and alpha sweep
    ├── analysis_window.m         ← Windowing function comparison
    └── generate_comparison_table.m ← Full system waveform × detector matrix
```

---

## How to Run

**Requirements:** MATLAB R2020a or later (no additional toolboxes needed)

### Interactive Simulation
```matlab
main_radar
```
The script will prompt you step by step:
- Waveform type: `LFM` / `PhaseCode` / `Hybrid`
- Window function: `rectangular` / `hamming` / `hanning` / `blackman`
- Clutter filters: MTI, Adaptive, Doppler (yes/no each)
- Detection method: `fixed` / `CA` / `OS` / `GO`
- Swerling model: `0` to `4`
- Filter type: `matched` / `mismatched_bw` / `mismatched_window`
- Number of targets + range, amplitude, velocity for each

### Analysis Scripts (run independently)
```matlab
plot_ambiguity           % 2D ambiguity function surfaces and slices
analysis_SNR_vs_Pd       % Monte Carlo SNR sweep
analysis_bandwidth       % Bandwidth vs resolution trade-off
analysis_CFAR_params     % CFAR parameter sensitivity
analysis_window          % Window function comparison
generate_comparison_table % Full waveform × detector benchmark
```

---

## System Parameters

| Parameter | Value |
|---|---|
| Sampling frequency | 10 MHz |
| Pulse duration | 10 µs |
| Bandwidth | 1 MHz (default) |
| Carrier frequency | 1 GHz |
| Speed of light | 3×10⁸ m/s |
| Range cell size | 15 m |
| Theoretical range resolution (c/2B) | 150 m |

---

## Key Features

### Waveforms
| Waveform | Time-Bandwidth Product | PSL | Range-Doppler Coupling |
|---|---|---|---|
| LFM Chirp | B×T = 10 | −13 dB | Yes (diagonal ridge) |
| Barker-13 Phase Code | 13 | −22 dB | No |
| Hybrid (LFM + Barker-13) | B×T×13 | −13 to −22 dB | Minimal |

### Detection Methods
- **Fixed threshold** — baseline, sensitive to noise level changes
- **CA-CFAR** — averages all training cells; optimal in uniform clutter
- **OS-CFAR** — uses kth-order statistic; robust to interfering targets
- **GO-CFAR** — takes greater of left/right averages; conservative in clutter edges

### Clutter Filters
- **MTI (3-pulse canceller)** — subtracts consecutive pulses to cancel stationary clutter
- **Adaptive LMS** — learns and subtracts clutter estimate (µ = 0.001)
- **Doppler notch** — zeros out ±1000 Hz around zero Doppler in frequency domain

### Target Models
| Swerling Model | Fluctuation | Distribution |
|---|---|---|
| 0 | None | Constant amplitude |
| 1 | Slow (per scan) | Rayleigh |
| 2 | Fast (per pulse) | Rayleigh |
| 3 | Slow (per scan) | Chi-squared (4 DOF) |
| 4 | Fast (per pulse) | Chi-squared (4 DOF) |

---

## Analysis Results Summary

- **Minimum SNR for Pd ≥ 0.90:** CA-CFAR achieves this earliest among adaptive detectors
- **Best range resolution:** LFM and Hybrid (theoretical c/2B = 150 m at 1 MHz)
- **Lowest sidelobes:** Blackman window (−58 dB PSL) at the cost of ~3 dB detection loss
- **Best two-target separation:** Hamming/Blackman windows suppress MF sidelobes enough to reveal weak targets 200 m behind a strong one
- **CFAR sweet spot:** alpha ≈ 2.0–2.5, training cells = 16 per side balances Pd and Pfa across all three CFAR variants

---

## Outputs

Each run of `main_radar.m` produces two figures:
1. **Waveform & Processing** — transmitted signal, received signal before/after filtering, matched filter output (full + zoomed), FFT spectrum, SNR before vs after
2. **Detection Results** — MF output with threshold overlay, detected target markers, Pd/Pfa bar chart

Each analysis script produces its own multi-panel figure saved with descriptive titles.

---

## Author

**Muhammad Asjad** | Roll No. 2023386
CE363 Digital Signal Processing — Complex Engineering Problem
Department of Computer Engineering, GIKI
