% =========================================================
%  generate_comparison_table.m
%  CE363 - Digital Signal Processing | Complex Engineering Problem
% =========================================================

clc; clear; close all;

fprintf('================================================\n');
fprintf('  FULL SYSTEM COMPARISON TABLE\n');
fprintf('================================================\n\n');

% ---------------------------------------------------------
%  FIXED PARAMETERS
% ---------------------------------------------------------
fs            = 10e6;
T             = 10e-6;
B             = 1e6;
c             = 3e8;
SNR_dB        = 10;
N_monte_carlo = 100;

num_train = 16;
num_guard = 2;
alpha     = 2.0;
k_os      = round(0.60 * 2 * num_train);   % BUG 1 FIX: 60th percentile = 19

tolerance_m   = 50;
target_ranges = [500, 800];
target_amps   = [0.8, 0.7];

waveforms = {'LFM', 'PhaseCode', 'Hybrid'};
detectors = {'fixed', 'CA', 'OS', 'GO'};
num_wav   = length(waveforms);
num_det   = length(detectors);

results_Pd        = zeros(num_wav, num_det);
results_Pfa       = zeros(num_wav, num_det);
results_snr_gain  = zeros(num_wav, num_det);
results_range_err = NaN(num_wav, num_det);
results_time      = zeros(num_wav, num_det);
psl_per_wav       = zeros(1, num_wav);

% ---------------------------------------------------------
%  PSL PER WAVEFORM
% ---------------------------------------------------------
for w = 1:num_wav
    switch waveforms{w}
        case 'LFM',       s_ref = generate_LFM(T, B, fs);
        case 'PhaseCode', s_ref = generate_PhaseCode(fs, T);
        case 'Hybrid',    s_ref = generate_Hybrid(T, B, fs);
    end
    N_p   = length(s_ref);
    ac    = abs(xcorr(s_ref, s_ref));
    ac    = ac(ceil(end/2):end);
    ac    = ac / max(ac + eps);
    ac_dB = 20*log10(ac + eps);
    mask  = max(10, round(N_p * 0.25));   % BUG 5 FIX: wider mask
    psl_per_wav(w) = compute_PSL_tbl(ac_dB, mask);
end

% ---------------------------------------------------------
%  MAIN LOOP
% ---------------------------------------------------------
fprintf('[Running %d configurations x %d trials...]\n\n', num_wav*num_det, N_monte_carlo);

cfg = 0;
for w = 1:num_wav
    wav_name = waveforms{w};

    switch wav_name
        case 'LFM',       transmitted = generate_LFM(T, B, fs);
        case 'PhaseCode', transmitted = generate_PhaseCode(fs, T);
        case 'Hybrid',    transmitted = generate_Hybrid(T, B, fs);
    end
    N_p = length(transmitted);

    % Analytical SNR gain
    switch wav_name
        case 'LFM',       tbp = B * T;
        case 'PhaseCode', tbp = 13;
        case 'Hybrid',    tbp = B * T * 13;
    end
    snr_gain_theory = 10 * log10(tbp);

    % Build clean echo
    max_d_s  = round(2 * max(target_ranges) * fs / c);
    recv_len = N_p + max_d_s + round(fs * 50e-6);
    echo_clean = zeros(1, recv_len);
    for tgt = 1:length(target_ranges)
        ds = round(2 * target_ranges(tgt) * fs / c);
        echo_clean(ds+1:ds+N_p) = echo_clean(ds+1:ds+N_p) + target_amps(tgt) * transmitted;
    end

    sp = mean(abs(echo_clean).^2);
    if sp == 0, sp = 1e-6; end
    np = sp / 10^(SNR_dB/10);

    for d = 1:num_det
        det_name = detectors{d};
        cfg = cfg + 1;
        fprintf('[%d/12] Waveform=%-10s Detection=%-8s ...', cfg, wav_name, det_name);

        t0      = tic;
        pd_acc  = 0;
        pfa_acc = 0;
        re_acc  = 0;
        re_cnt  = 0;

        for trial = 1:N_monte_carlo
            noise    = sqrt(np/2) * (randn(1,recv_len) + 1j*randn(1,recv_len));
            received = echo_clean + noise;

            mf       = abs(xcorr(received, transmitted));
            mf       = mf(ceil(end/2):end);
            range_ax = (0:length(mf)-1) * c / (2*fs);
            N_cells  = length(mf);

            win    = hamming(N_cells)';
            mf_win = mf .* win;

            % BUG 1: k_os=19, BUG 2: GO alpha=1.5, BUG 3: fixed uses median*4, BUG 4: CA alpha=2.5
            switch det_name
                case 'fixed'
                    noise_est = median(mf_win) * 4.0;   % BUG 3 FIX
                    [di,~] = detect_fixed(mf_win, noise_est);
                case 'CA'
                    [di,~] = detect_CA_CFAR(mf_win, num_train, num_guard, 2.5);   % BUG 4 FIX
                case 'OS'
                    [di,~] = detect_OS_CFAR(mf_win, num_train, num_guard, alpha, k_os);
                case 'GO'
                    [di,~] = detect_GO_CFAR(mf_win, num_train, num_guard, 1.5);   % BUG 2 FIX
            end

            dr = get_ranges_tbl(di, range_ax, mf_win);
            [pd, pfa, td, ~] = score_tbl(target_ranges, dr, N_cells, tolerance_m);

            pd_acc  = pd_acc  + pd;
            pfa_acc = pfa_acc + pfa;

            if td > 0
                for tgt = 1:length(target_ranges)
                    diffs = abs(dr - target_ranges(tgt));
                    [mn, ~] = min(diffs);
                    if mn <= tolerance_m
                        re_acc = re_acc + mn;
                        re_cnt = re_cnt + 1;
                    end
                end
            end
        end

        results_Pd(w,d)       = pd_acc  / N_monte_carlo;
        results_Pfa(w,d)      = pfa_acc / N_monte_carlo;
        results_snr_gain(w,d) = snr_gain_theory;

        if re_cnt > 0
            results_range_err(w,d) = re_acc / re_cnt;
        end

        results_time(w,d) = toc(t0);

        fprintf('  Pd=%.2f Pfa=%.4f SNR_gain=%.1fdB (%.1fs)\n', ...
            results_Pd(w,d), results_Pfa(w,d), results_snr_gain(w,d), results_time(w,d));
    end
end

% ---------------------------------------------------------
%  PRINT TABLE
% ---------------------------------------------------------
fprintf('\n============================================================\n');
fprintf('  FULL SYSTEM COMPARISON TABLE\n');
fprintf('============================================================\n');
fprintf('%-12s %-10s %-12s %-8s %-14s %-12s\n', ...
    'Waveform','Detector','SNR_gain(dB)','Pd','Pfa','RangeErr(m)');
fprintf('%s\n', repmat('-',1,72));
for w = 1:num_wav
    for d = 1:num_det
        if isnan(results_range_err(w,d))
            rstr = 'N/A';
        else
            rstr = sprintf('%.1f', results_range_err(w,d));
        end
        fprintf('%-12s %-10s %-12.1f %-8.2f %-14.5f %-12s\n', ...
            waveforms{w}, detectors{d}, results_snr_gain(w,d), ...
            results_Pd(w,d), results_Pfa(w,d), rstr);
    end
end
fprintf('============================================================\n\n');

% ---------------------------------------------------------
%  BEST CONFIGURATION SUMMARY
% ---------------------------------------------------------
fprintf('--- Best Configuration Per Metric ---\n');

valid = results_Pd >= 0.5;

[bv,bi] = max(results_Pd(:));
[bw,bd] = ind2sub([num_wav num_det], bi);
fprintf('  Highest Pd       : %s + %s = %.2f\n', waveforms{bw}, detectors{bd}, bv);

tp = results_Pfa; tp(~valid) = inf;
[bv,bi] = min(tp(:));
if bv < inf
    [bw,bd] = ind2sub([num_wav num_det], bi);
    fprintf('  Lowest Pfa       : %s + %s = %.5f\n', waveforms{bw}, detectors{bd}, bv);
else
    fprintf('  Lowest Pfa       : No valid config (Pd>=0.5 required)\n');
end

[bv,bi] = max(results_snr_gain(:));
[bw,bd] = ind2sub([num_wav num_det], bi);
fprintf('  Highest SNR gain : %s + %s = %.1f dB\n', waveforms{bw}, detectors{bd}, bv);

te = results_range_err; te(~valid | isnan(results_range_err)) = inf;
[bv,bi] = min(te(:));
if bv < inf
    [bw,bd] = ind2sub([num_wav num_det], bi);
    fprintf('  Best range accur.: %s + %s = %.1f m\n', waveforms{bw}, detectors{bd}, bv);
else
    fprintf('  Best range accur.: No valid config (Pd>=0.5 required)\n');
end
fprintf('\n');

% ---------------------------------------------------------
%  PLOTS
% ---------------------------------------------------------
figure('Name','System Comparison','NumberTitle','off','Position',[50 50 1400 900]);

det_cats = categorical(detectors, detectors);

subplot(2,3,1);
bar(det_cats, results_Pd', 0.8);
yline(0.9,'--k','LineWidth',1,'HandleVisibility','off');
title('Pd — All Configurations');
xlabel('Detector'); ylabel('Pd'); ylim([0 1.15]); grid on;
legend(waveforms,'Location','southwest');

subplot(2,3,2);
bar(det_cats, results_Pfa', 0.8);
title('Pfa — All Configurations');
xlabel('Detector'); ylabel('Pfa'); grid on;
legend(waveforms,'Location','northeast');

subplot(2,3,3);
bar(det_cats, results_snr_gain', 0.8);
title('SNR Gain (dB) — All Configurations');
xlabel('Detector'); ylabel('SNR Gain (dB)'); grid on;
legend(waveforms,'Location','northeast');

subplot(2,3,4);
bar(categorical(waveforms, waveforms), psl_per_wav, 0.5, 'FaceColor',[0.7 0.5 0.3]);
yline(-13,'--b','LineWidth',1,'HandleVisibility','off');
title('Peak Sidelobe Level per Waveform');
xlabel('Waveform'); ylabel('PSL (dB)'); grid on;

subplot(2,3,5);
re_plot = results_range_err;
re_plot(isnan(re_plot)) = 0;
bar(det_cats, re_plot', 0.8);
title('Range Error (0 = N/A)');
xlabel('Detector'); ylabel('Range Error (m)'); grid on;
legend(waveforms,'Location','northeast');

subplot(2,3,6);
overall_score = results_Pd - 10*results_Pfa;
overall_score(overall_score < 0) = 0;
imagesc(overall_score);
axis xy;
xticks(1:num_det); xticklabels(detectors);
yticks(1:num_wav);  yticklabels(waveforms);
colorbar;
title('Overall Score (Pd - 10xPfa)');
xlabel('Detector'); ylabel('Waveform');
for w2 = 1:num_wav
    for d2 = 1:num_det
        text(d2, w2, sprintf('%.2f',overall_score(w2,d2)), ...
            'HorizontalAlignment','center','FontSize',9,'Color','w','FontWeight','bold');
    end
end

% =========================================================
%  HELPER FUNCTIONS
% =========================================================

function psl_dB = compute_PSL_tbl(slice_dB, mask_hw)
    [~,pi] = max(slice_dB);
    N = length(slice_dB);
    masked = slice_dB;
    masked(max(1,pi-mask_hw):min(N,pi+mask_hw)) = -999;
    ps = max(masked);
    if ps <= -999, psl_dB = -60; else, psl_dB = ps; end
end

function dr = get_ranges_tbl(indices, range_axis, mf)
    dr = [];
    if isempty(indices), return; end
    clusters = {}; cluster = indices(1);
    for i = 2:length(indices)
        if indices(i)-indices(i-1) <= 10
            cluster = [cluster, indices(i)];
        else
            clusters{end+1} = cluster; cluster = indices(i);
        end
    end
    clusters{end+1} = cluster;
    for i = 1:length(clusters)
        cv = clusters{i};
        [~,pk] = max(mf(cv));
        dr(end+1) = range_axis(cv(pk));
    end
end

function [Pd, Pfa, td, fd] = score_tbl(true_r, det_r, total_cells, tol)
    nt = numel(true_r); nd = numel(det_r);
    tm = false(1,nt); dm = false(1,nd);
    for i = 1:nt
        ix = find(abs(det_r-true_r(i))<=tol & ~dm, 1,'first');
        if ~isempty(ix), tm(i)=true; dm(ix)=true; end
    end
    td  = sum(tm); fd = nd-sum(dm);
    Pd  = td/max(nt,1);
    Pfa = fd/max(total_cells-nt*20,1);
end