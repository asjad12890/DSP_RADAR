% =========================================================
%  plot_ambiguity.m
%  CE363 - Digital Signal Processing | Complex Engineering Problem
% =========================================================

clc; clear; close all;

% ---------------------------------------------------------
%  PARAMETERS
% ---------------------------------------------------------
fs          = 10e6;
T           = 10e-6;
B           = 1e6;
c           = 3e8;
num_doppler = 129;      % odd so centre = exact 0 Hz
max_doppler = 200000;
pulse_len   = round(T * fs);   % 100 samples
dx_m        = c / (2 * fs);    % 15 m per sample
dB_floor    = -40;

fprintf('Waveforms generated.\n');
fprintf('Pulse length: %d samples\n', pulse_len);
fprintf('Range cell size: %.1f m\n', dx_m);
fprintf('Theoretical range resolution: %.1f m\n', c/(2*B));
fprintf('Theoretical Doppler resolution: %.1f Hz\n', 1/T);

% ---------------------------------------------------------
%  GENERATE WAVEFORMS
% ---------------------------------------------------------
s_LFM       = generate_LFM(T, B, fs);
s_PhaseCode = generate_PhaseCode(fs, T);
s_Hybr_raw  = generate_Hybrid(T, B, fs);

% Trim to pulse_len
s_LFM       = s_LFM(1:pulse_len);
s_PhaseCode = s_PhaseCode(1:pulse_len);
s_Hybr_raw  = s_Hybr_raw(1:pulse_len);

% BUG 2 FIX: check if generate_Hybrid is identical to PhaseCode
similarity = abs(dot(s_Hybr_raw/norm(s_Hybr_raw), s_PhaseCode/norm(s_PhaseCode)));
if similarity > 0.99
    fprintf('WARNING: generate_Hybrid identical to PhaseCode — using manual Hybrid construction\n');
    barker13      = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1];
    chip_len      = floor(pulse_len / 13);
    phase_code_up = repelem(barker13, chip_len);
    phase_code_up = phase_code_up(1:pulse_len);
    s_Hybr        = s_LFM .* phase_code_up;
    s_Hybr        = s_Hybr / norm(s_Hybr);
else
    s_Hybr = s_Hybr_raw;
end

% ---------------------------------------------------------
%  COMPUTE 2D AMBIGUITY FUNCTIONS
% ---------------------------------------------------------
delay_axis   = (0:pulse_len-1) * dx_m;
doppler_axis = linspace(-max_doppler, max_doppler, num_doppler);
t_vec        = (0:pulse_len-1) / fs;

AF_LFM       = zeros(num_doppler, pulse_len);
AF_PhaseCode = zeros(num_doppler, pulse_len);
AF_Hybrid    = zeros(num_doppler, pulse_len);

fprintf('Computing ambiguity functions...\n');

for d = 1:num_doppler
    dp = exp(1j * 2*pi * doppler_axis(d) * t_vec);

    corr_out = xcorr(s_LFM.*dp, s_LFM);
    corr_out = corr_out(ceil(end/2):end);
    AF_LFM(d,:) = abs(corr_out(1:pulse_len));

    corr_out = xcorr(s_PhaseCode.*dp, s_PhaseCode);
    corr_out = corr_out(ceil(end/2):end);
    AF_PhaseCode(d,:) = abs(corr_out(1:pulse_len));

    corr_out = xcorr(s_Hybr.*dp, s_Hybr);
    corr_out = corr_out(ceil(end/2):end);
    AF_Hybrid(d,:) = abs(corr_out(1:pulse_len));
end

% Normalise by global peak
AF_LFM       = AF_LFM       / max(max(AF_LFM));
AF_PhaseCode = AF_PhaseCode / max(max(AF_PhaseCode));
AF_Hybrid    = AF_Hybrid    / max(max(AF_Hybrid));

% Convert to dB
AF_LFM_dB       = 20*log10(AF_LFM       + eps);
AF_PhaseCode_dB = 20*log10(AF_PhaseCode + eps);
AF_Hybrid_dB    = 20*log10(AF_Hybrid    + eps);

fprintf('Ambiguity functions computed.\n\n');

% ---------------------------------------------------------
%  EXTRACT SLICES
% ---------------------------------------------------------
zero_row = ceil(num_doppler/2);   % exact fd=0 row

% Zero-Doppler slices (range profile at fd=0)
slice_LFM       = AF_LFM_dB(zero_row, :);
slice_PhaseCode = AF_PhaseCode_dB(zero_row, :);
slice_Hybrid    = AF_Hybrid_dB(zero_row, :);

% Zero-delay slices (Doppler profile at tau=0)
dslice_LFM       = AF_LFM_dB(:,1)';
dslice_PhaseCode = AF_PhaseCode_dB(:,1)';
dslice_Hybrid    = AF_Hybrid_dB(:,1)';
doppler_axis_plot = doppler_axis;

% ---------------------------------------------------------
%  BUG 1 FIX: MEASURE ALL METRICS FROM ACTUAL DATA
% ---------------------------------------------------------

% --- Range resolution from zero-Doppler slice ---
res_LFM       = measure_resolution_m(10.^(slice_LFM/20),       delay_axis);
res_PhaseCode = measure_resolution_m(10.^(slice_PhaseCode/20), delay_axis);
res_Hybrid    = measure_resolution_m(10.^(slice_Hybrid/20),    delay_axis);

% --- Doppler resolution from zero-delay slice ---
dres_LFM       = measure_resolution_hz(10.^(dslice_LFM/20),       doppler_axis_plot);
dres_PhaseCode = measure_resolution_hz(10.^(dslice_PhaseCode/20), doppler_axis_plot);
dres_Hybrid    = measure_resolution_hz(10.^(dslice_Hybrid/20),    doppler_axis_plot);

% --- PSL from zero-Doppler slice ---
% mask based on resolution: mask_hw = round(res/dx_m) + 3
mask_LFM   = max(4, round(res_LFM/dx_m)       + 3);
mask_PC    = max(4, round(res_PhaseCode/dx_m)  + 3);
mask_Hybr  = max(4, round(res_Hybrid/dx_m)     + 3);

psl_LFM       = compute_PSL(slice_LFM,       mask_LFM);
psl_PhaseCode = compute_PSL(slice_PhaseCode, mask_PC);
psl_Hybrid    = compute_PSL(slice_Hybrid,    mask_Hybr);

% --- Range-Doppler coupling: measure ridge tilt angle ---
coupling_LFM       = measure_coupling(AF_LFM,       delay_axis, doppler_axis);
coupling_PhaseCode = measure_coupling(AF_PhaseCode, delay_axis, doppler_axis);
coupling_Hybrid    = measure_coupling(AF_Hybrid,    delay_axis, doppler_axis);

thresh_slope = 500;   % Hz/m threshold for "coupled"
if abs(coupling_LFM) > thresh_slope
    coupling_str_LFM = sprintf('YES (slope=%.0fk Hz/m)', coupling_LFM/1000);
else
    coupling_str_LFM = sprintf('NO  (slope=%.0f Hz/m)',  coupling_LFM);
end
if abs(coupling_PhaseCode) > thresh_slope
    coupling_str_PC = sprintf('YES (slope=%.0fk Hz/m)', coupling_PhaseCode/1000);
else
    coupling_str_PC = sprintf('NO  (slope=%.0f Hz/m)',  coupling_PhaseCode);
end
if abs(coupling_Hybrid) > thresh_slope
    coupling_str_Hybr = sprintf('YES (slope=%.0fk Hz/m)', coupling_Hybrid/1000);
else
    coupling_str_Hybr = sprintf('NO  (slope=%.0f Hz/m)',  coupling_Hybrid);
end

% --- Print metrics table ---
fprintf('============================================\n');
fprintf('  AMBIGUITY FUNCTION METRICS\n');
fprintf('============================================\n');
fprintf('%-28s %-14s %-14s %-14s\n', 'Metric','LFM','PhaseCode','Hybrid');
fprintf('%-28s %-14.1f %-14.1f %-14.1f\n', 'Range res. -3dB (m)',       res_LFM,       res_PhaseCode, res_Hybrid);
fprintf('%-28s %-14.1f %-14.1f %-14.1f\n', 'Doppler res. -3dB (Hz)',    dres_LFM,      dres_PhaseCode,dres_Hybrid);
fprintf('%-28s %-14.1f %-14.1f %-14.1f\n', 'Peak Sidelobe Level (dB)',  psl_LFM,       psl_PhaseCode, psl_Hybrid);
fprintf('%-28s %-14s %-14s %-14s\n',        'Range-Doppler coupling',    coupling_str_LFM, coupling_str_PC, coupling_str_Hybr);
fprintf('============================================\n\n');

% ---------------------------------------------------------
%  BUG 3 FIX: mirror zero-Doppler slices for symmetric display
% ---------------------------------------------------------
slice_LFM_full       = [fliplr(slice_LFM(2:end)),       slice_LFM];
slice_PhaseCode_full = [fliplr(slice_PhaseCode(2:end)), slice_PhaseCode];
slice_Hybrid_full    = [fliplr(slice_Hybrid(2:end)),    slice_Hybrid];
delay_axis_full      = [-fliplr(delay_axis(2:end)),     delay_axis];

% ---------------------------------------------------------
%  FIGURE 1 — 3D AMBIGUITY SURFACES
% ---------------------------------------------------------
figure('Name','3D Ambiguity Functions','NumberTitle','off','Position',[50 50 1600 500]);

[D_mesh, F_mesh] = meshgrid(delay_axis/1000, doppler_axis);

subplot(1,3,1);
AF_plot = max(AF_LFM_dB, dB_floor);
surf(D_mesh, F_mesh, AF_plot, 'EdgeColor','none'); shading interp;
clim([dB_floor 0]);
colorbar; colormap(jet);
xlabel('Delay / Range (km)'); ylabel('Doppler (Hz)'); zlabel('Magnitude (dB)');
title('LFM Chirp — Ambiguity Function');
view(20, 55);   % BUG 5: LFM uses special angle to show diagonal ridge

subplot(1,3,2);
AF_plot = max(AF_PhaseCode_dB, dB_floor);
surf(D_mesh, F_mesh, AF_plot, 'EdgeColor','none'); shading interp;
clim([dB_floor 0]);
colorbar; colormap(jet);
xlabel('Delay / Range (km)'); ylabel('Doppler (Hz)'); zlabel('Magnitude (dB)');
title('Phase-Coded (Barker-13) — Ambiguity Function');
view(45, 35);

subplot(1,3,3);
AF_plot = max(AF_Hybrid_dB, dB_floor);
surf(D_mesh, F_mesh, AF_plot, 'EdgeColor','none'); shading interp;
clim([dB_floor 0]);
colorbar; colormap(jet);
xlabel('Delay / Range (km)'); ylabel('Doppler (Hz)'); zlabel('Magnitude (dB)');
title('Hybrid — Ambiguity Function');
view(45, 35);

sgtitle('DSP CEP, M.ASJAD, 2023386')


% ---------------------------------------------------------
%  FIGURE 2 — TOP-DOWN CONTOUR VIEW
% ---------------------------------------------------------
figure('Name','Top-Down View','NumberTitle','off','Position',[50 600 1600 500]);

subplot(1,3,1);
imagesc(delay_axis/1000, doppler_axis, max(AF_LFM_dB, dB_floor));
set(gca,'YDir','normal'); colorbar; colormap(jet); clim([dB_floor 0]);
xlabel('Delay / Range (km)'); ylabel('Doppler (Hz)');
title('LFM — Top View');
xlim([0 0.3]);   % BUG 6: zoom to 300m

subplot(1,3,2);
imagesc(delay_axis/1000, doppler_axis, max(AF_PhaseCode_dB, dB_floor));
set(gca,'YDir','normal'); colorbar; colormap(jet); clim([dB_floor 0]);
xlabel('Delay / Range (km)'); ylabel('Doppler (Hz)');
title('Phase-Coded — Top View');
xlim([0 0.3]);

subplot(1,3,3);
imagesc(delay_axis/1000, doppler_axis, max(AF_Hybrid_dB, dB_floor));
set(gca,'YDir','normal'); colorbar; colormap(jet); clim([dB_floor 0]);
xlabel('Delay / Range (km)'); ylabel('Doppler (Hz)');
title('Hybrid — Top View');
xlim([0 0.3]);

sgtitle('DSP CEP, M.ASJAD, 2023386')


% ---------------------------------------------------------
%  FIGURE 3 — ZERO-DOPPLER SLICE (BUG 3+7 FIX)
% ---------------------------------------------------------
figure('Name','Zero-Doppler Slice','NumberTitle','off','Position',[50 50 900 500]);

% BUG 7: correct plot order so legend matches
plot(delay_axis_full/1000, slice_LFM_full,       'b-',  'LineWidth', 2.0); hold on;
plot(delay_axis_full/1000, slice_PhaseCode_full, 'r-',  'LineWidth', 2.5);
plot(delay_axis_full/1000, slice_Hybrid_full,    'g--', 'LineWidth', 2.0);
yline(-3,  '--k', '-3 dB', 'LineWidth', 1);
yline(-13, '--m', '-13 dB (rect sidelobe)', 'LineWidth', 1);
title('Zero-Doppler Slice — Range Resolution and Sidelobe Comparison');
xlabel('Range Delay (m)'); ylabel('Normalised Amplitude (dB)');
legend('LFM', 'PhaseCode (Barker-13)', 'Hybrid', 'Location', 'northeast');
grid on; ylim([-50 5]);
xlim([-500 500]/1000);   % BUG 7: zoom to mainlobe region

sgtitle('DSP CEP, M.ASJAD, 2023386')


% ---------------------------------------------------------
%  FIGURE 4 — ZERO-DELAY SLICE (BUG 4 FIX)
% ---------------------------------------------------------
figure('Name','Zero-Delay Slice','NumberTitle','off','Position',[50 550 900 500]);

% BUG 4: plot PhaseCode first (bottom), then Hybrid, then LFM on top
plot(doppler_axis_plot, dslice_PhaseCode, 'r',   'LineWidth', 2.5); hold on;
plot(doppler_axis_plot, dslice_Hybrid,    'g--', 'LineWidth', 2.0);
plot(doppler_axis_plot, dslice_LFM,       'b',   'LineWidth', 2.0);
yline(-3, '--k', '-3dB', 'LineWidth', 1);
title('Zero-Delay Slice — Doppler Tolerance Comparison');
xlabel('Doppler Frequency (Hz)'); ylabel('Normalised Amplitude (dB)');
legend('PhaseCode (Barker-13)', 'Hybrid', 'LFM', 'Location', 'northeast');
grid on; ylim([-50 5]);
xlim([-150000 150000]);   % BUG 4: zoom to ±1.5× Doppler resolution

sgtitle('DSP CEP, M.ASJAD, 2023386')


% =========================================================
%  HELPER FUNCTIONS
% =========================================================

function width_m = measure_resolution_m(ac_linear, axis_m)
    % Measure -3dB width from one-sided autocorrelation in metres
    % ac_linear normalised so peak=1, axis_m is delay in metres
    threshold = 0.5;
    [~, peak_idx] = max(ac_linear);
    N = length(ac_linear);

    right_idx = peak_idx;
    while right_idx < N && ac_linear(right_idx) >= threshold
        right_idx = right_idx + 1;
    end

    if right_idx > 1 && right_idx <= N
        v_hi  = ac_linear(right_idx-1);
        v_lo  = ac_linear(right_idx);
        frac  = (v_hi - threshold) / max(v_hi - v_lo, eps);
        right_cross_m = axis_m(right_idx-1) + frac*(axis_m(min(right_idx,N)) - axis_m(right_idx-1));
    else
        right_cross_m = axis_m(min(right_idx, N));
    end

    if peak_idx <= 2
        width_m = 2 * (right_cross_m - axis_m(1));
    else
        left_idx = peak_idx;
        while left_idx > 1 && ac_linear(left_idx) >= threshold
            left_idx = left_idx - 1;
        end
        if left_idx >= 1 && left_idx < peak_idx
            v_hi  = ac_linear(left_idx+1);
            v_lo  = ac_linear(left_idx);
            frac  = (v_hi - threshold) / max(v_hi - v_lo, eps);
            left_cross_m = axis_m(left_idx+1) - frac*(axis_m(left_idx+1) - axis_m(max(left_idx,1)));
        else
            left_cross_m = axis_m(left_idx);
        end
        width_m = right_cross_m - left_cross_m;
    end
    width_m = max(width_m, axis_m(2) - axis_m(1));
end

function width_hz = measure_resolution_hz(ac_linear, axis_hz)
    % Measure -3dB width from Doppler slice in Hz
    % Doppler slice is symmetric, peak at centre
    threshold = 0.5;
    [~, peak_idx] = max(ac_linear);
    N = length(ac_linear);

    right_idx = peak_idx;
    while right_idx < N && ac_linear(right_idx) >= threshold
        right_idx = right_idx + 1;
    end
    if right_idx > 1 && right_idx <= N
        v_hi = ac_linear(right_idx-1); v_lo = ac_linear(right_idx);
        frac = (v_hi - threshold) / max(v_hi - v_lo, eps);
        right_hz = axis_hz(right_idx-1) + frac*(axis_hz(min(right_idx,N)) - axis_hz(right_idx-1));
    else
        right_hz = axis_hz(min(right_idx,N));
    end

    left_idx = peak_idx;
    while left_idx > 1 && ac_linear(left_idx) >= threshold
        left_idx = left_idx - 1;
    end
    if left_idx >= 1 && left_idx < peak_idx
        v_hi = ac_linear(left_idx+1); v_lo = ac_linear(left_idx);
        frac = (v_hi - threshold) / max(v_hi - v_lo, eps);
        left_hz = axis_hz(left_idx+1) - frac*(axis_hz(left_idx+1) - axis_hz(max(left_idx,1)));
    else
        left_hz = axis_hz(left_idx);
    end

    width_hz = max(right_hz - left_hz, abs(axis_hz(2) - axis_hz(1)));
end

function psl_dB = compute_PSL(slice_dB, mask_hw)
    % Measure PSL from dB slice with mask_hw samples blanked each side of peak
    [~, peak_idx] = max(slice_dB);
    N = length(slice_dB);
    masked = slice_dB;
    masked(max(1,peak_idx-mask_hw):min(N,peak_idx+mask_hw)) = -999;
    peak_sl = max(masked);
    if peak_sl <= -999
        psl_dB = -60;
    else
        psl_dB = peak_sl;
    end
end

function slope = measure_coupling(AF, delay_axis, doppler_axis)
    % Measure ridge tilt: for each Doppler row find peak delay location
    % fit line through (delay_peak, doppler) pairs
    % return slope in Hz/m
    [num_dop, ~] = size(AF);
    peak_delays = zeros(1, num_dop);
    for d = 1:num_dop
        [~, idx] = max(AF(d,:));
        peak_delays(d) = delay_axis(idx);
    end
    % polyfit: doppler = slope * delay + offset
    p = polyfit(peak_delays, doppler_axis, 1);
    slope = p(1);   % Hz/m
end