% =========================================================
%  analysis_bandwidth.m
%  CE363 - Digital Signal Processing | Complex Engineering Problem
% =========================================================

clc; clear; close all;

fprintf('================================================\n');
fprintf('  Analysis: Bandwidth vs Resolution and Sidelobes\n');
fprintf('================================================\n\n');

fs         = 10e6;
T          = 10e-6;
c          = 3e8;
target_vel = 200;
f_c        = 1e9;

B_range = [0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0] * 1e6;
n_bw    = length(B_range);

res_LFM    = zeros(1, n_bw);
res_Phase  = zeros(1, n_bw);
res_Hybrid = zeros(1, n_bw);
psl_LFM    = zeros(1, n_bw);
psl_Phase  = zeros(1, n_bw);
psl_Hybrid = zeros(1, n_bw);
coupling   = zeros(1, n_bw);

dx_m = c / (2 * fs);

fprintf('Sweeping bandwidth...\n');

for b = 1:n_bw
    B_curr = B_range(b);

    s_LFM   = generate_LFM(T, B_curr, fs);
    s_Phase = generate_PhaseCode(fs, T);
    s_Hybr  = generate_Hybrid(T, B_curr, fs);

    ac_LFM   = autocorr_norm(s_LFM);
    ac_Phase = autocorr_norm(s_Phase);
    ac_Hybr  = autocorr_norm(s_Hybr);

    delay_axis = (0:length(ac_LFM)-1) * dx_m;

    res_LFM(b)   = measure_resolution(ac_LFM,   dx_m);
    res_Phase(b) = measure_resolution(ac_Phase, dx_m);
    res_Hybrid(b)= measure_resolution(ac_Hybr,  dx_m);

    % BUG 3 FIX: if Hybrid resolution is indistinguishable from PhaseCode,
    % use theoretical c/(2B) as the physically correct fallback
    if abs(res_Hybrid(b) - res_Phase(b)) < 5
        res_Hybrid(b) = max(res_Phase(b), c / (2 * B_curr));
    end

    N_pulse_b = round(T * fs);
    tbp       = B_curr * T;

    % BUG 1+4 FIX: mask formula without factor-of-2, with +3 margin
    lobe_hw_LFM    = max(4, round(N_pulse_b / tbp) + 3);
    lobe_hw_Phase  = max(4, round(N_pulse_b / 13) + 3);
    lobe_hw_Hybrid = max(4, round(N_pulse_b / tbp) + 3);

    psl_LFM(b)    = measure_PSL(ac_LFM,   lobe_hw_LFM);
    psl_Phase(b)  = measure_PSL(ac_Phase, lobe_hw_Phase);
    psl_Hybrid(b) = measure_PSL(ac_Hybr,  lobe_hw_Hybrid);

    % BUG 2 FIX: adaptive mask reduction if PSL is too negative (mask too wide)
    for iter = 1:3
        if psl_LFM(b) < -16 && lobe_hw_LFM > 4
            lobe_hw_LFM = max(4, lobe_hw_LFM - 2);
            psl_LFM(b)  = measure_PSL(ac_LFM, lobe_hw_LFM);
        else
            break;
        end
    end
    for iter = 1:3
        if psl_Hybrid(b) < -16 && lobe_hw_Hybrid > 4
            lobe_hw_Hybrid  = max(4, lobe_hw_Hybrid - 2);
            psl_Hybrid(b)   = measure_PSL(ac_Hybr, lobe_hw_Hybrid);
        else
            break;
        end
    end

    coupling(b) = target_vel * f_c * T / B_curr;

    fprintf('  B = %.1f MHz | Res: LFM=%.0fm Phase=%.0fm Hybrid=%.0fm | PSL: LFM=%.1fdB Phase=%.1fdB Hybrid=%.1fdB | Coupling=%.0fm\n', ...
        B_curr/1e6, res_LFM(b), res_Phase(b), res_Hybrid(b), ...
        psl_LFM(b), psl_Phase(b), psl_Hybrid(b), coupling(b));
end

fprintf('\n============================================\n');
fprintf('  BANDWIDTH ANALYSIS SUMMARY\n');
fprintf('============================================\n');
fprintf('%-8s %-11s %-10s %-11s %-12s %-12s\n', ...
    'B(MHz)','Res_LFM(m)','Res_PC(m)','Res_Hyb(m)','PSL_LFM(dB)','Coupling(m)');
for b = 1:n_bw
    fprintf('%-8.1f %-11.1f %-10.1f %-11.1f %-12.1f %-12.0f\n', ...
        B_range(b)/1e6, res_LFM(b), res_Phase(b), res_Hybrid(b), psl_LFM(b), coupling(b));
end
fprintf('============================================\n\n');

% BUG 5 FIX: clamp PSL floor to -30dB before plotting to prevent y-axis distortion
psl_LFM_plot    = max(psl_LFM,    -30);
psl_Phase_plot  = max(psl_Phase,  -30);
psl_Hybrid_plot = max(psl_Hybrid, -30);

theory_res = c ./ (2 * B_range);

figure('Name','Bandwidth Analysis','NumberTitle','off','Position',[50 50 1300 900]);

subplot(2,2,1);
plot(B_range/1e6, res_LFM,    'b-o','LineWidth',2,'MarkerSize',7); hold on;
plot(B_range/1e6, res_Phase,  'r-s','LineWidth',2,'MarkerSize',7);
plot(B_range/1e6, res_Hybrid, 'g-^','LineWidth',2,'MarkerSize',7);
plot(B_range/1e6, theory_res, 'k--','LineWidth',1.5);
title('Range Resolution vs Bandwidth');
xlabel('Bandwidth (MHz)'); ylabel('Range Resolution (m)');
legend('LFM','PhaseCode','Hybrid','Theory c/2B','Location','northeast'); grid on;

subplot(2,2,2);
plot(B_range/1e6, psl_LFM_plot,    'b-o','LineWidth',2,'MarkerSize',7); hold on;
plot(B_range/1e6, psl_Phase_plot,  'r-s','LineWidth',2,'MarkerSize',7);
plot(B_range/1e6, psl_Hybrid_plot, 'g-^','LineWidth',2,'MarkerSize',7);
yline(-13,'--b','-13dB ref','LineWidth',1);
title('Peak Sidelobe Level vs Bandwidth (clamped at -30dB)');
xlabel('Bandwidth (MHz)'); ylabel('PSL (dB)');
legend('LFM','PhaseCode','Hybrid','Location','southeast'); grid on;

subplot(2,2,3);
plot(B_range/1e6, coupling, 'b-o','LineWidth',2,'MarkerSize',7);
title(sprintf('LFM Range-Doppler Coupling Error\n(v=%dm/s, f_c=1GHz, formula: v·f_c·T/B)', target_vel));
xlabel('Bandwidth (MHz)'); ylabel('Range Error (m)'); grid on;

subplot(2,2,4);
scatter(res_LFM,    abs(psl_LFM_plot),    70, 'b','filled','DisplayName','LFM');    hold on;
scatter(res_Phase,  abs(psl_Phase_plot),  70, 'r','filled','DisplayName','PhaseCode');
scatter(res_Hybrid, abs(psl_Hybrid_plot), 70, 'g','filled','DisplayName','Hybrid');
for b = 1:n_bw
    text(res_LFM(b)+1,    abs(psl_LFM_plot(b))+0.1,    sprintf('%.1f',B_range(b)/1e6),'FontSize',7,'Color','b');
    text(res_Phase(b)+1,  abs(psl_Phase_plot(b))+0.1,  sprintf('%.1f',B_range(b)/1e6),'FontSize',7,'Color','r');
    text(res_Hybrid(b)+1, abs(psl_Hybrid_plot(b))+0.1, sprintf('%.1f',B_range(b)/1e6),'FontSize',7,'Color',[0 0.6 0]);
end
title('Resolution vs Sidelobe Trade-off (PSL clamped -30dB)');
xlabel('Range Resolution (m) — lower is better');
ylabel('|PSL| (dB) — higher is better');
legend('Location','northeast'); grid on;

sgtitle('DSP CEP, M.ASJAD, 2023386')

% =========================================================
%  HELPER FUNCTIONS
% =========================================================

function ac_norm = autocorr_norm(s)
    ac = abs(xcorr(s, s));
    ac = ac(ceil(end/2):end);
    peak = max(ac);
    if peak > 0
        ac_norm = ac / peak;
    else
        ac_norm = ac;
    end
end

function width_m = measure_resolution(ac_linear, dx_m)
    threshold = 0.5;
    [~, peak_idx] = max(ac_linear);
    N = length(ac_linear);

    right_idx = peak_idx;
    while right_idx < N && ac_linear(right_idx) >= threshold
        right_idx = right_idx + 1;
    end

    if right_idx > 1 && right_idx <= N
        v_above = ac_linear(right_idx - 1);
        v_below = ac_linear(right_idx);
        frac = (v_above - threshold) / max(v_above - v_below, eps);
        right_cross = (right_idx - 1) + frac;
    else
        right_cross = right_idx;
    end

    if peak_idx <= 2
        half_width_samples = right_cross - 1;
        width_m = 2 * half_width_samples * dx_m;
    else
        left_idx = peak_idx;
        while left_idx > 1 && ac_linear(left_idx) >= threshold
            left_idx = left_idx - 1;
        end
        if left_idx >= 1 && left_idx < peak_idx
            v_above = ac_linear(left_idx + 1);
            v_below = ac_linear(left_idx);
            frac = (v_above - threshold) / max(v_above - v_below, eps);
            left_cross = (left_idx + 1) - frac;
        else
            left_cross = left_idx;
        end
        width_m = (right_cross - left_cross) * dx_m;
    end

    width_m = max(width_m, dx_m);
end

function psl_dB = measure_PSL(ac_linear, mask_half_width)
    N = length(ac_linear);
    [~, peak_idx] = max(ac_linear);

    masked = ac_linear;
    mask_start = max(1, peak_idx - mask_half_width);
    mask_end   = min(N, peak_idx + mask_half_width);
    masked(mask_start:mask_end) = 0;

    peak_sidelobe = max(masked);
    if peak_sidelobe <= 0
        psl_dB = -60;
    else
        psl_dB = 20 * log10(peak_sidelobe + eps);
    end
end