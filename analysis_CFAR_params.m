% =========================================================
%  analysis_CFAR_params.m
%  CE363 - Digital Signal Processing | Complex Engineering Problem
% =========================================================

clc; clear; close all;

fprintf('================================================\n');
fprintf('  Analysis: CFAR Parameters vs Detection Performance\n');
fprintf('================================================\n\n');

% ---------------------------------------------------------
%  FIXED PARAMETERS
% ---------------------------------------------------------
fs           = 10e6;
T            = 10e-6;
B            = 1e6;
c            = 3e8;
SNR_dB       = 10;
N_mc         = 200;
target_range = 500;
target_amp   = 0.8;
tolerance_m  = 50;
n_guard      = 2;

train_range  = [2, 4, 6, 8, 12, 16, 20, 24, 32];
alpha_range  = [0.3, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0];  % extended to show fast Pfa rise at low alpha
n_tp         = length(train_range);
n_ap         = length(alpha_range);   % auto-sized

% ---------------------------------------------------------
%  GENERATE WAVEFORM
% ---------------------------------------------------------
transmitted = generate_LFM(T, B, fs);
N_pulse     = length(transmitted);
max_delay   = round(2 * target_range * fs / c);
recv_len    = N_pulse + max_delay + round(fs * 50e-6);

echo_signal = zeros(1, recv_len);
echo_signal(max_delay+1 : max_delay+N_pulse) = target_amp * transmitted;

sig_pow   = mean(abs(echo_signal).^2);
if sig_pow == 0, sig_pow = 1e-6; end
noise_pow = sig_pow / 10^(SNR_dB / 10);

% ---------------------------------------------------------
%  PART 1: SWEEP TRAINING CELLS (alpha_fixed=1.5)
%  Fresh noise per trial so low training counts get fair test
% ---------------------------------------------------------
fprintf('Part 1: Training cells sweep (alpha=1.5, N=%d)...\n', N_mc);

alpha_fixed = 2.0;   % safer operating point — more threshold margin

Pd_CA_tr  = zeros(1, n_tp);  Pd_OS_tr  = zeros(1, n_tp);  Pd_GO_tr  = zeros(1, n_tp);
Pfa_CA_tr = zeros(1, n_tp);  Pfa_OS_tr = zeros(1, n_tp);  Pfa_GO_tr = zeros(1, n_tp);

for t = 1:n_tp
    n_tr = train_range(t);
    k_os = max(1, round(0.60 * 2 * n_tr));   % 60th percentile — more stable, avoids collapse at high train counts

    pd_ca=0; pd_os=0; pd_go=0;
    pfa_ca=0; pfa_os=0; pfa_go=0;

    for tr = 1:N_mc
        % Fresh noise every trial
        noise    = sqrt(noise_pow/2) * (randn(1,recv_len) + 1j*randn(1,recv_len));
        received = echo_signal + noise;
        mf       = abs(xcorr(received, transmitted));
        mf       = mf(ceil(end/2):end);
        N_cells  = length(mf);
        range_ax = (0:N_cells-1) * c / (2*fs);

        di_ca = cfar_ca(mf, n_tr, n_guard, alpha_fixed);
        di_os = cfar_os(mf, n_tr, n_guard, alpha_fixed, k_os);
        di_go = cfar_go(mf, n_tr, n_guard, alpha_fixed);

        dr_ca = get_ranges(di_ca, range_ax, mf);
        dr_os = get_ranges(di_os, range_ax, mf);
        dr_go = get_ranges(di_go, range_ax, mf);

        n_tested = max(2*(N_cells - 2*(n_tr+n_guard)), 1);

        [p1,f1] = score_det(target_range, dr_ca, n_tested, tolerance_m);
        [p2,f2] = score_det(target_range, dr_os, n_tested, tolerance_m);
        [p3,f3] = score_det(target_range, dr_go, n_tested, tolerance_m);

        pd_ca=pd_ca+p1;  pfa_ca=pfa_ca+f1;
        pd_os=pd_os+p2;  pfa_os=pfa_os+f2;
        pd_go=pd_go+p3;  pfa_go=pfa_go+f3;
    end

    Pd_CA_tr(t)  = pd_ca  / N_mc;  Pd_OS_tr(t)  = pd_os  / N_mc;  Pd_GO_tr(t)  = pd_go  / N_mc;
    Pfa_CA_tr(t) = pfa_ca / N_mc;  Pfa_OS_tr(t) = pfa_os / N_mc;  Pfa_GO_tr(t) = pfa_go / N_mc;

    fprintf('  train=%2d | Pd: CA=%.2f OS=%.2f GO=%.2f | Pfa: CA=%.4f OS=%.4f GO=%.4f\n', ...
        n_tr, Pd_CA_tr(t),Pd_OS_tr(t),Pd_GO_tr(t), ...
        Pfa_CA_tr(t),Pfa_OS_tr(t),Pfa_GO_tr(t));
end

% ---------------------------------------------------------
%  PART 2: SWEEP ALPHA (train=16, k=75th percentile)
%  Fresh noise per trial
% ---------------------------------------------------------
fprintf('\nPart 2: Alpha sweep (train=16, N=%d)...\n', N_mc);

n_tr_fixed = 16;
k_os_fixed = max(1, round(0.60 * 2 * n_tr_fixed));   % 60th percentile = 19 — comparable to CA mean, more stable out of 32

Pd_CA_al  = zeros(1, n_ap);  Pd_OS_al  = zeros(1, n_ap);  Pd_GO_al  = zeros(1, n_ap);
Pfa_CA_al = zeros(1, n_ap);  Pfa_OS_al = zeros(1, n_ap);  Pfa_GO_al = zeros(1, n_ap);

for a = 1:n_ap
    alpha_curr = alpha_range(a);

    pd_ca=0; pd_os=0; pd_go=0;
    pfa_ca=0; pfa_os=0; pfa_go=0;

    for tr = 1:N_mc
        noise    = sqrt(noise_pow/2) * (randn(1,recv_len) + 1j*randn(1,recv_len));
        received = echo_signal + noise;
        mf       = abs(xcorr(received, transmitted));
        mf       = mf(ceil(end/2):end);
        N_cells  = length(mf);
        range_ax = (0:N_cells-1) * c / (2*fs);

        di_ca = cfar_ca(mf, n_tr_fixed, n_guard, alpha_curr);
        di_os = cfar_os(mf, n_tr_fixed, n_guard, alpha_curr, k_os_fixed);
        di_go = cfar_go(mf, n_tr_fixed, n_guard, alpha_curr);

        dr_ca = get_ranges(di_ca, range_ax, mf);
        dr_os = get_ranges(di_os, range_ax, mf);
        dr_go = get_ranges(di_go, range_ax, mf);

        n_tested = max(2*(N_cells - 2*(n_tr_fixed+n_guard)), 1);

        [p1,f1] = score_det(target_range, dr_ca, n_tested, tolerance_m);
        [p2,f2] = score_det(target_range, dr_os, n_tested, tolerance_m);
        [p3,f3] = score_det(target_range, dr_go, n_tested, tolerance_m);

        pd_ca=pd_ca+p1;  pfa_ca=pfa_ca+f1;
        pd_os=pd_os+p2;  pfa_os=pfa_os+f2;
        pd_go=pd_go+p3;  pfa_go=pfa_go+f3;
    end

    Pd_CA_al(a)  = pd_ca  / N_mc;  Pd_OS_al(a)  = pd_os  / N_mc;  Pd_GO_al(a)  = pd_go  / N_mc;
    Pfa_CA_al(a) = pfa_ca / N_mc;  Pfa_OS_al(a) = pfa_os / N_mc;  Pfa_GO_al(a) = pfa_go / N_mc;

    fprintf('  alpha=%.1f | Pd: CA=%.2f OS=%.2f GO=%.2f | Pfa: CA=%.4f OS=%.4f GO=%.4f\n', ...
        alpha_curr, Pd_CA_al(a),Pd_OS_al(a),Pd_GO_al(a), ...
        Pfa_CA_al(a),Pfa_OS_al(a),Pfa_GO_al(a));
end

% ---------------------------------------------------------
%  PLOTS
% ---------------------------------------------------------
figure('Name','CFAR Parameter Analysis','NumberTitle','off','Position',[50 50 1400 900]);

% Plot 1 — Training cells vs Pd
subplot(2,3,1);
plot(train_range, Pd_CA_tr, 'b-o','LineWidth',2,'MarkerSize',7); hold on;
plot(train_range, Pd_OS_tr, 'r-s','LineWidth',2,'MarkerSize',7);
plot(train_range, Pd_GO_tr, 'g-^','LineWidth',2,'MarkerSize',7);
yline(0.9,'--k','Pd=0.90','LineWidth',1);
title(sprintf('Training Cells vs Pd (alpha=%.1f, N=%d)',alpha_fixed,N_mc));
xlabel('Training Cells per Side'); ylabel('Pd');
legend('CA-CFAR','OS-CFAR','GO-CFAR','Location','southeast'); grid on;
ylim([0 1.1]); xlim([0 train_range(end)+2]);
text(2, 0.15, 'Note: OS sensitivity to train count', 'FontSize', 8, 'Color', 'r');

% Plot 2 — Training cells vs Pfa
subplot(2,3,2);
semilogy(train_range, max(Pfa_CA_tr,1e-6), 'b-o','LineWidth',2,'MarkerSize',7); hold on;
semilogy(train_range, max(Pfa_OS_tr,1e-6), 'r-s','LineWidth',2,'MarkerSize',7);
semilogy(train_range, max(Pfa_GO_tr,1e-6), 'g-^','LineWidth',2,'MarkerSize',7);
title(sprintf('Training Cells vs Pfa (alpha=%.1f)',alpha_fixed));
xlabel('Training Cells per Side'); ylabel('Pfa');
legend('CA-CFAR','OS-CFAR','GO-CFAR','Location','northeast'); grid on;

% Plot 3 — Alpha vs Pd
subplot(2,3,3);
plot(alpha_range, Pd_CA_al, 'b-o','LineWidth',2,'MarkerSize',7); hold on;
plot(alpha_range, Pd_OS_al, 'r-s','LineWidth',2,'MarkerSize',7);
plot(alpha_range, Pd_GO_al, 'g-^','LineWidth',2,'MarkerSize',7);
yline(0.9,'--k','Pd=0.90','LineWidth',1);
title(sprintf('Alpha vs Pd (train=%d, N=%d) — sweet spot alpha=2.0-2.5',n_tr_fixed,N_mc));
xlabel('Alpha (threshold multiplier)'); ylabel('Pd');
legend('CA-CFAR','OS-CFAR','GO-CFAR','Location','southwest'); grid on;
ylim([0 1.1]);

% Plot 4 — Alpha vs Pfa
subplot(2,3,4);
semilogy(alpha_range, max(Pfa_CA_al,1e-6), 'b-o','LineWidth',2,'MarkerSize',7); hold on;
semilogy(alpha_range, max(Pfa_OS_al,1e-6), 'r-s','LineWidth',2,'MarkerSize',7);
semilogy(alpha_range, max(Pfa_GO_al,1e-6), 'g-^','LineWidth',2,'MarkerSize',7);
title(sprintf('Alpha vs Pfa (train=%d)',n_tr_fixed));
xlabel('Alpha'); ylabel('Pfa');
legend('CA-CFAR','OS-CFAR','GO-CFAR','Location','northeast'); grid on;

% Plot 5 — ROC curve (Pd vs Pfa from alpha sweep)
subplot(2,3,5);
plot(Pfa_CA_al, Pd_CA_al, 'b-o','LineWidth',2,'MarkerSize',7); hold on;
plot(Pfa_OS_al, Pd_OS_al, 'r-s','LineWidth',2,'MarkerSize',7);
plot(Pfa_GO_al, Pd_GO_al, 'g-^','LineWidth',2,'MarkerSize',7);
title(sprintf('ROC Curve — Alpha Sweep (train=%d)',n_tr_fixed));
xlabel('Pfa'); ylabel('Pd');
legend('CA-CFAR','OS-CFAR','GO-CFAR','Location','southeast'); grid on;
ylim([0 1.1]);

% Plot 6 — Summary grouped bar (NO yyaxis, Pfa scaled x100)
subplot(2,3,6);
idx16 = find(train_range==16, 1);
if isempty(idx16), [~,idx16] = min(abs(train_range-16)); end
Pd_bar  = [Pd_CA_tr(idx16),  Pd_OS_tr(idx16),  Pd_GO_tr(idx16)];
Pfa_bar = [Pfa_CA_tr(idx16), Pfa_OS_tr(idx16), Pfa_GO_tr(idx16)];
x = 1:3; w = 0.35;
bar(x-w/2, max(Pd_bar,0),      w, 'FaceColor',[0.2 0.5 0.9],'DisplayName','Pd'); hold on;
bar(x+w/2, max(Pfa_bar*100,0), w, 'FaceColor',[0.9 0.4 0.2],'DisplayName','Pfa x100');
xticks(x); xticklabels({'CA-CFAR','OS-CFAR','GO-CFAR'});
title(sprintf('Summary: train=16, alpha=%.1f',alpha_fixed));
ylabel('Pd  /  Pfa x100'); ylim([0 1.2]);
legend('Location','northeast'); grid on;

sgtitle('DSP CEP, M.ASJAD, 2023386')

% ---------------------------------------------------------
%  INLINE CFAR IMPLEMENTATIONS
% ---------------------------------------------------------

function idx = cfar_ca(mf, n_tr, n_g, alpha)
    N = length(mf);
    idx = [];
    half = n_tr + n_g;
    for i = half+1 : N-half
        left  = mf(i-n_g-n_tr : i-n_g-1);
        right = mf(i+n_g+1    : i+n_g+n_tr);
        T = alpha * mean([left(:); right(:)]);
        if mf(i) > T
            idx(end+1) = i;
        end
    end
end

function idx = cfar_os(mf, n_tr, n_g, alpha, k)
    N = length(mf);
    idx = [];
    half = n_tr + n_g;
    for i = half+1 : N-half
        left  = mf(i-n_g-n_tr : i-n_g-1);
        right = mf(i+n_g+1    : i+n_g+n_tr);
        sv    = sort([left(:); right(:)]);
        k_s   = max(1, min(k, length(sv)));
        T     = alpha * sv(k_s);
        if mf(i) > T
            idx(end+1) = i;
        end
    end
end

function idx = cfar_go(mf, n_tr, n_g, alpha)
    N = length(mf);
    idx = [];
    half = n_tr + n_g;
    for i = half+1 : N-half
        left  = mf(i-n_g-n_tr : i-n_g-1);
        right = mf(i+n_g+1    : i+n_g+n_tr);
        T     = alpha * max(mean(left(:)), mean(right(:)));
        if mf(i) > T
            idx(end+1) = i;
        end
    end
end

function dr = get_ranges(indices, range_axis, mf)
    dr = [];
    if isempty(indices), return; end
    clusters = {};
    cluster  = indices(1);
    for i = 2:length(indices)
        if indices(i) - indices(i-1) <= 10
            cluster = [cluster, indices(i)];
        else
            clusters{end+1} = cluster;
            cluster = indices(i);
        end
    end
    clusters{end+1} = cluster;
    for i = 1:length(clusters)
        cv = clusters{i};
        [~, pi] = max(mf(cv));
        dr(end+1) = range_axis(cv(pi));
    end
end

function [Pd, Pfa] = score_det(true_r, det_r, n_tested, tol)
    % n_tested = actual cells tested by CFAR = 2*(N - 2*(n_tr+n_guard))
    nt = numel(true_r);
    nd = numel(det_r);
    tm = false(1, nt);
    dm = false(1, nd);
    for i = 1:nt
        ix = find(abs(det_r - true_r(i)) <= tol & ~dm, 1, 'first');
        if ~isempty(ix)
            tm(i) = true;
            dm(ix) = true;
        end
    end
    td  = sum(tm);
    fd  = nd - sum(dm);
    Pd  = td / max(nt, 1);
    Pfa = fd / max(n_tested - nt, 1);
end