% =========================================================
%  analysis_SNR_vs_Pd.m
%  CE363 - Digital Signal Processing | Complex Engineering Problem
%
%  Parametric Analysis: SNR vs Probability of Detection
%
%  What this script does:
%  - Sweeps SNR from -10 dB to 20 dB
%  - At each SNR runs multiple Monte Carlo trials
%    (repeat with different random noise each time and average)
%  - Measures Pd and Pfa for all 4 detection methods
%  - Plots detection curves for comparison
%  - Prints minimum SNR needed for Pd >= 0.9 for each method
% =========================================================

clc; clear; close all;

fprintf('================================================\n');
fprintf('  Analysis: SNR vs Probability of Detection\n');
fprintf('================================================\n\n');

% ---------------------------------------------------------
%  FIXED PARAMETERS
% ---------------------------------------------------------

fs            = 10e6;     % sampling frequency
T             = 10e-6;    % pulse duration
B             = 1e6;      % bandwidth
c             = 3e8;      % speed of light
fc            = 1e9;      % carrier frequency
N_monte_carlo = 50;       % number of trials per SNR point
                          % higher = more accurate but slower
                          % 50 gives good accuracy in reasonable time

% single target at 500m, amplitude 0.8, stationary
target_range  = 500;
target_amp    = 0.8;
target_vel    = 0;

% CFAR parameters
num_train     = 8;
num_guard     = 2;
alpha         = 3.0;
k_os          = 12;
fixed_factor  = 3.0;
tolerance_m   = 50;       % detection within 50m of true range = correct

% SNR range to sweep
snr_range     = -10 : 2 : 20;   % from -10 to 20 dB in steps of 2
num_snr       = length(snr_range);

% ---------------------------------------------------------
%  PRE-ALLOCATE RESULTS ARRAYS
% ---------------------------------------------------------

Pd_fixed = zeros(1, num_snr);
Pd_CA    = zeros(1, num_snr);
Pd_OS    = zeros(1, num_snr);
Pd_GO    = zeros(1, num_snr);

Pfa_fixed = zeros(1, num_snr);
Pfa_CA    = zeros(1, num_snr);
Pfa_OS    = zeros(1, num_snr);
Pfa_GO    = zeros(1, num_snr);

% ---------------------------------------------------------
%  GENERATE TRANSMITTED WAVEFORM (LFM — done once)
% ---------------------------------------------------------

transmitted  = generate_LFM(T, B, fs);
pulse_length = length(transmitted);

% ---------------------------------------------------------
%  MAIN LOOP — SWEEP SNR
% ---------------------------------------------------------

fprintf('Running Monte Carlo analysis...\n');
fprintf('SNR range: %d to %d dB | Trials per point: %d\n\n', ...
        snr_range(1), snr_range(end), N_monte_carlo);

for s = 1 : num_snr

    SNR_dB = snr_range(s);

    % accumulators for this SNR level
    det_fixed = 0;  fa_fixed = 0;
    det_CA    = 0;  fa_CA    = 0;
    det_OS    = 0;  fa_OS    = 0;
    det_GO    = 0;  fa_GO    = 0;

    for trial = 1 : N_monte_carlo

        % --- build received signal ---
        max_delay     = round(2 * target_range * fs / c);
        recv_len      = pulse_length + max_delay + round(fs * 50e-6);
        received      = zeros(1, recv_len);

        % place echo
        delay_samp    = round(2 * target_range * fs / c);
        received(delay_samp+1 : delay_samp+pulse_length) = target_amp * transmitted;

        % add AWGN noise
        sig_pow       = mean(abs(received).^2);
        if sig_pow == 0, sig_pow = 1e-6; end
        noise_pow     = sig_pow / 10^(SNR_dB/10);
        noise         = sqrt(noise_pow/2) * (randn(1,recv_len) + 1j*randn(1,recv_len));
        received      = received + noise;

        % add simple clutter
        t_r           = (0:recv_len-1)/fs;
        clutter       = 0.05 * exp(1j*2*pi*0.04e6*t_r) .* (0.3 + 0.2*randn(1,recv_len));
        received      = received + clutter;

        % --- matched filter ---
        mf            = abs(xcorr(received, transmitted));
        N_c           = length(mf);
        mf            = mf(ceil(N_c/2):end);
        range_ax      = (0:length(mf)-1) * c / (2*fs);

        % --- run all 4 detection methods ---
        [di_f, ~] = detect_fixed(mf,   fixed_factor);
        [di_ca,~] = detect_CA_CFAR(mf, num_train, num_guard, alpha);
        [di_os,~] = detect_OS_CFAR(mf, num_train, num_guard, alpha, k_os);
        [di_go,~] = detect_GO_CFAR(mf, num_train, num_guard, alpha);

        % convert detected indices to ranges
        dr_f  = get_detected_ranges(di_f,  range_ax, mf);
        dr_ca = get_detected_ranges(di_ca, range_ax, mf);
        dr_os = get_detected_ranges(di_os, range_ax, mf);
        dr_go = get_detected_ranges(di_go, range_ax, mf);

        % score each method
        [pd_f, pfa_f, td_f, fd_f]   = compute_detection_metrics([target_range], dr_f,  length(mf), tolerance_m);
        [pd_ca,pfa_ca,td_ca,fd_ca]  = compute_detection_metrics([target_range], dr_ca, length(mf), tolerance_m);
        [pd_os,pfa_os,td_os,fd_os]  = compute_detection_metrics([target_range], dr_os, length(mf), tolerance_m);
        [pd_go,pfa_go,td_go,fd_go]  = compute_detection_metrics([target_range], dr_go, length(mf), tolerance_m);

        det_fixed = det_fixed + pd_f;   fa_fixed = fa_fixed + pfa_f;
        det_CA    = det_CA    + pd_ca;  fa_CA    = fa_CA    + pfa_ca;
        det_OS    = det_OS    + pd_os;  fa_OS    = fa_OS    + pfa_os;
        det_GO    = det_GO    + pd_go;  fa_GO    = fa_GO    + pfa_go;

    end

    % average over all trials
    Pd_fixed(s)  = det_fixed / N_monte_carlo;
    Pd_CA(s)     = det_CA    / N_monte_carlo;
    Pd_OS(s)     = det_OS    / N_monte_carlo;
    Pd_GO(s)     = det_GO    / N_monte_carlo;

    Pfa_fixed(s) = fa_fixed  / N_monte_carlo;
    Pfa_CA(s)    = fa_CA     / N_monte_carlo;
    Pfa_OS(s)    = fa_OS     / N_monte_carlo;
    Pfa_GO(s)    = fa_GO     / N_monte_carlo;

    fprintf('  SNR = %+4d dB | Pd: Fixed=%.2f CA=%.2f OS=%.2f GO=%.2f\n', ...
            SNR_dB, Pd_fixed(s), Pd_CA(s), Pd_OS(s), Pd_GO(s));
end

% ---------------------------------------------------------
%  FIND MINIMUM SNR FOR Pd >= 0.9
% ---------------------------------------------------------

fprintf('\n--- Minimum SNR for Pd >= 0.90 ---\n');
methods  = {'Fixed', 'CA-CFAR', 'OS-CFAR', 'GO-CFAR'};
pd_all   = {Pd_fixed, Pd_CA, Pd_OS, Pd_GO};

for m = 1:4
    idx = find(pd_all{m} >= 0.9, 1, 'first');
    if isempty(idx)
        fprintf('  %-10s : never reaches Pd=0.90 in tested range\n', methods{m});
    else
        fprintf('  %-10s : Pd >= 0.90 at SNR = %+d dB\n', methods{m}, snr_range(idx));
    end
end

% ---------------------------------------------------------
%  PLOTS
% ---------------------------------------------------------

figure('Name','SNR vs Pd — Detection Method Comparison', ...
       'NumberTitle','off','Position',[100 100 1100 500]);

% Plot 1 — Pd vs SNR
subplot(1,2,1);
plot(snr_range, Pd_fixed, 'r-o',  'LineWidth',1.5,'MarkerSize',5); hold on;
plot(snr_range, Pd_CA,    'b-s',  'LineWidth',1.5,'MarkerSize',5);
plot(snr_range, Pd_OS,    'g-^',  'LineWidth',1.5,'MarkerSize',5);
plot(snr_range, Pd_GO,    'm-d',  'LineWidth',1.5,'MarkerSize',5);
yline(0.9, '--k', 'Pd = 0.90 target', 'LineWidth', 1.2);
title('Probability of Detection vs SNR');
xlabel('SNR (dB)'); ylabel('Pd');
legend('Fixed','CA-CFAR','OS-CFAR','GO-CFAR','Location','southeast');
ylim([0 1.05]); grid on;

% Plot 2 — Pfa vs SNR
subplot(1,2,2);
semilogy(snr_range, Pfa_fixed+eps, 'r-o', 'LineWidth',1.5,'MarkerSize',5); hold on;
semilogy(snr_range, Pfa_CA+eps,    'b-s', 'LineWidth',1.5,'MarkerSize',5);
semilogy(snr_range, Pfa_OS+eps,    'g-^', 'LineWidth',1.5,'MarkerSize',5);
semilogy(snr_range, Pfa_GO+eps,    'm-d', 'LineWidth',1.5,'MarkerSize',5);
title('Probability of False Alarm vs SNR');
xlabel('SNR (dB)'); ylabel('Pfa (log scale)');
legend('Fixed','CA-CFAR','OS-CFAR','GO-CFAR','Location','northeast');
grid on;

% ---------------------------------------------------------
%  HELPER FUNCTIONS
% ---------------------------------------------------------

function detected_ranges = get_detected_ranges(detected_indices, range_axis, mf_output)
% clusters nearby detections and returns one range per cluster
    detected_ranges = [];
    if isempty(detected_indices), return; end

    clusters = {};
    cluster  = detected_indices(1);
    for i = 2:length(detected_indices)
        if detected_indices(i) - detected_indices(i-1) <= 3
            cluster = [cluster, detected_indices(i)];
        else
            clusters{end+1} = cluster;
            cluster = detected_indices(i);
        end
    end
    clusters{end+1} = cluster;

    for i = 1:length(clusters)
        c = clusters{i};
        [~, peak_idx] = max(mf_output(c));
        detected_ranges(end+1) = range_axis(c(peak_idx));
    end
end

function [Pd, Pfa, true_detected, false_detections] = compute_detection_metrics( ...
         true_ranges, detected_ranges, total_cells, tolerance_m)
    num_tgt           = numel(true_ranges);
    num_det           = numel(detected_ranges);
    target_matched    = false(1, num_tgt);
    detection_matched = false(1, num_det);

    for i = 1:num_tgt
        idx = find(abs(detected_ranges - true_ranges(i)) <= tolerance_m ...
                   & ~detection_matched, 1, 'first');
        if ~isempty(idx)
            target_matched(i)    = true;
            detection_matched(idx) = true;
        end
    end

    true_detected    = sum(target_matched);
    false_detections = num_det - sum(detection_matched);
    Pd               = true_detected / max(num_tgt, 1);
    Pfa              = false_detections / max(total_cells - num_tgt, 1);
end
