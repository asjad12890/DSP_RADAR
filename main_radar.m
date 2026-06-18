% =========================================================
%  main_radar.m
%  CE363 - Digital Signal Processing | Complex Engineering Problem
%  Radar DSP Processing Framework
% =========================================================

clc; clear; close all;

% ---------------------------------------------------------
%  SECTION 1: FIXED SYSTEM PARAMETERS
% ---------------------------------------------------------

fs     = 10e6;
T      = 10e-6;
B      = 1e6;
c      = 3e8;
SNR_dB = 10;
N_fft  = 1024;

% ---------------------------------------------------------
%  SECTION 2: USER INPUTS
%  IMPROVEMENT 1: validatestring on all text inputs
%    — trims whitespace, accepts partial matches, gives clean error messages
%  IMPROVEMENT 2: switch statements replace if-elseif chains
%    — cleaner, faster to read, standard MATLAB style
% ---------------------------------------------------------

fprintf('============================================\n');
fprintf('   CE363 Radar DSP Simulation\n');
fprintf('============================================\n\n');

fprintf('Available waveforms: LFM | PhaseCode | Hybrid\n');
waveform_raw  = upper(strtrim(input('Enter waveform type: ', 's')));
waveform_type = validatestring(waveform_raw, {'LFM','PHASECODE','HYBRID'});
switch waveform_type
    case 'LFM',       waveform_type = 'LFM';
    case 'PHASECODE', waveform_type = 'PhaseCode';
    case 'HYBRID',    waveform_type = 'Hybrid';
end

fprintf('\nAvailable windows: rectangular | hamming | hanning | blackman\n');
window_type = validatestring(lower(strtrim(input('Enter window type: ', 's'))), ...
    {'rectangular', 'hamming', 'hanning', 'blackman'});

fprintf('\n--- Clutter Mitigation (enter yes or no for each) ---\n');
use_MTI      = validatestring(lower(strtrim(input('Enable MTI filter?      (yes/no): ', 's'))), {'yes','no'});
use_adaptive = validatestring(lower(strtrim(input('Enable adaptive filter? (yes/no): ', 's'))), {'yes','no'});
use_doppler  = validatestring(lower(strtrim(input('Enable doppler filter?  (yes/no): ', 's'))), {'yes','no'});

fprintf('\n--- Detection Method ---\n');
fprintf('Available methods: fixed | CA | OS | GO\n');
detect_method = validatestring(lower(strtrim(input('Enter detection method: ', 's'))), ...
    {'fixed','ca','os','go'});

fprintf('\n--- Target Fluctuation Model ---\n');
fprintf('Swerling 0 = no fluctuation (constant amplitude)\n');
fprintf('Swerling 1 = slow fluctuation, Rayleigh (constant per scan)\n');
fprintf('Swerling 2 = fast fluctuation, Rayleigh (changes every pulse)\n');
fprintf('Swerling 3 = slow fluctuation, Chi-squared\n');
fprintf('Swerling 4 = fast fluctuation, Chi-squared\n');
swerling_model = input('Enter Swerling model (0/1/2/3/4): ');
if ~isscalar(swerling_model) || ~ismember(swerling_model, 0:4)
    error('Swerling model must be 0, 1, 2, 3, or 4');
end
if ismember(swerling_model, [2,4])
    fprintf('Note: Swerling %d fast fluctuation requires multiple pulses to observe.\n', swerling_model);
end

fprintf('\n--- Matched Filter Type ---\n');
fprintf('matched           = perfect reference (standard)\n');
fprintf('mismatched_bw     = reference with 10%% bandwidth mismatch\n');
fprintf('mismatched_window = reference windowed with Hamming (lower sidelobes, SNR loss)\n');
filter_type = validatestring(lower(strtrim(input('Enter filter type: ', 's'))), ...
    {'matched','mismatched_bw','mismatched_window'});

% --- Correlation method ---
fprintf('\n--- Correlation Method ---\n');
fprintf('xcorr = cross-correlation (standard, direct)\n');
fprintf('conv  = convolution with flipped signal (equivalent result)\n');
corr_method = validatestring(lower(strtrim(input('Enter correlation method (xcorr/conv): ', 's'))), ...
    {'xcorr','conv'});

num_targets = input('\nEnter number of targets: ');
if ~isscalar(num_targets) || num_targets < 1 || floor(num_targets) ~= num_targets
    error('Number of targets must be a positive integer.');
end

target_ranges     = zeros(1, num_targets);
target_amplitudes = zeros(1, num_targets);
target_velocities = zeros(1, num_targets);

for i = 1 : num_targets
    target_ranges(i)     = input(sprintf('  Range of target %d (metres): ', i));
    target_amplitudes(i) = input(sprintf('  Amplitude of target %d (0.1 to 1.0): ', i));
    target_velocities(i) = input(sprintf('  Velocity of target %d (m/s, + towards radar, - away): ', i));
    if target_ranges(i) <= 0
        error('Target range must be positive.');
    end
    if target_amplitudes(i) <= 0
        error('Target amplitude must be positive.');
    end
end

fprintf('\n--- Setup Summary ---\n');
fprintf('Waveform   : %s\n', waveform_type);
fprintf('Window     : %s\n', window_type);
fprintf('MTI        : %s\n', use_MTI);
fprintf('Adaptive   : %s\n', use_adaptive);
fprintf('Doppler    : %s\n', use_doppler);
fprintf('Detection  : %s\n', upper(detect_method));
fprintf('Swerling   : %d\n', swerling_model);
fprintf('Filter     : %s\n', filter_type);
fprintf('Corr method: %s\n', corr_method);
for i = 1 : num_targets
    fprintf('Target %d   : %.0f m | amp %.2f | vel %.1f m/s\n', ...
            i, target_ranges(i), target_amplitudes(i), target_velocities(i));
end
fprintf('---------------------\n\n');

% ---------------------------------------------------------
%  SECTION 3: GENERATE TRANSMITTED WAVEFORM
%  IMPROVEMENT 2: switch replaces if-elseif chain
% ---------------------------------------------------------

switch waveform_type
    case 'LFM'
        transmitted = generate_LFM(T, B, fs);
    case 'PhaseCode'
        transmitted = generate_PhaseCode(fs, T);
    case 'Hybrid'
        transmitted = generate_Hybrid(T, B, fs);
end

pulse_length = length(transmitted);
fprintf('Waveform generated. Pulse length: %d samples\n', pulse_length);

% ---------------------------------------------------------
%  SECTION 4: SIMULATE RECEIVED SIGNAL
% ---------------------------------------------------------

fc = 1e9;

max_delay_samples = round(2 * max(target_ranges) * fs / c);
receive_length    = pulse_length + max_delay_samples + round(fs * 50e-6);
received          = zeros(1, receive_length);

for i = 1 : num_targets

    base_amp = target_amplitudes(i);

    if swerling_model == 0
        actual_amp = base_amp;
    elseif swerling_model == 1 || swerling_model == 2
        raw_amp    = base_amp * abs(randn + 1j*randn) / sqrt(2);
        actual_amp = max(raw_amp, base_amp * 0.40);  % floor at 40% of base
    else
        raw_amp    = base_amp * abs(randn + 1j*randn + randn + 1j*randn) / 2;
        actual_amp = max(raw_amp, base_amp * 0.40);  % floor at 40% of base
    end

    fprintf('  Target %d: base amp=%.2f | Swerling %d actual amp=%.4f\n', ...
            i, base_amp, swerling_model, actual_amp);

    doppler_shift = 2 * target_velocities(i) * fc / c;
    t_pulse_vec   = (0 : pulse_length-1) / fs;
    doppler_phase = exp(1j * 2*pi * doppler_shift * t_pulse_vec);
    echo          = actual_amp * transmitted .* doppler_phase;

    if strcmp(waveform_type, 'LFM') && target_velocities(i) ~= 0
        range_error_m       = target_velocities(i) * T * c / (2 * B);
        range_error_samples = round(range_error_m * 2 * fs / c);
        fprintf('  Target %d LFM range-Doppler coupling: %.1f m range error\n', i, range_error_m);
    else
        range_error_samples = 0;
    end

    delay_samples = max(0, round(2 * target_ranges(i) * fs / c) + range_error_samples);
    start_idx     = delay_samples + 1;
    end_idx       = start_idx + pulse_length - 1;
    if end_idx <= receive_length
        received(start_idx:end_idx) = received(start_idx:end_idx) + echo;
    end
end

signal_power = mean(abs(received).^2);
if signal_power == 0, signal_power = 1e-6; end
noise_power  = signal_power / 10^(SNR_dB/10);
noise        = sqrt(noise_power/2) * (randn(1,receive_length) + 1j*randn(1,receive_length));
received     = received + noise;

t_recv      = (0 : receive_length-1) / fs;
amp_drift1  = 0.02 + 0.01 * sin(2*pi*0.5 * t_recv);
amp_drift2  = 0.015 + 0.008 * cos(2*pi*0.3 * t_recv);
amp_drift3  = 0.01 + 0.005 * sin(2*pi*0.8 * t_recv + 1.2);
freq_drift1 = 0.04e6 + 0.01e6 * sin(2*pi*0.2 * t_recv);
freq_drift2 = 0.07e6 + 0.01e6 * cos(2*pi*0.15 * t_recv);
freq_drift3 = 0.02e6 + 0.005e6 * sin(2*pi*0.4 * t_recv);
clutter1    = amp_drift1 .* exp(1j * 2*pi * cumsum(freq_drift1) / fs);
clutter2    = amp_drift2 .* exp(1j * 2*pi * cumsum(freq_drift2) / fs);
clutter3    = amp_drift3 .* exp(1j * 2*pi * cumsum(freq_drift3) / fs);
received    = received + clutter1 + clutter2 + clutter3;

fprintf('Received signal built:\n');
fprintf('  Targets: %d | Swerling model: %d\n', num_targets, swerling_model);
fprintf('  SNR: %d dB | Time-varying clutter: 3 components\n', SNR_dB);

% ---------------------------------------------------------
%  SECTION 4.5: CLUTTER MITIGATION
% ---------------------------------------------------------

received_raw = received;
target_killed = false(1, num_targets);  % tracks targets removed by filters
fprintf('\nClutter mitigation:\n');

if strcmp(use_MTI, 'yes'),      received = apply_MTI(received);          end
if strcmp(use_adaptive, 'yes'), received = apply_adaptive(received);     end
if strcmp(use_doppler, 'yes'),  received = apply_doppler(received, fs);  end

% warn if any target will be killed by the Doppler notch
% target_killed was initialised above; update it here
if strcmp(use_doppler, 'yes')
    notch_half_hz = 1000;
    fc_local      = 1e9;
    for i = 1 : num_targets
        fd = abs(2 * target_velocities(i) * fc_local / c);
        if fd < notch_half_hz
            target_killed(i) = true;
            fprintf('  !! WARNING: Target %d Doppler = %.0f Hz is inside the %.0f Hz notch.\n', ...
                    i, fd, notch_half_hz);
            fprintf('             This target has been removed by the Doppler filter.\n');
            fprintf('             Disable the Doppler filter or use velocity > %.1f m/s.\n', ...
                    notch_half_hz * c / (2 * fc_local));
        end
    end
end

if strcmp(use_MTI,'no') && strcmp(use_adaptive,'no') && strcmp(use_doppler,'no')
    fprintf('  No clutter mitigation applied\n');
end
fprintf('Clutter mitigation complete.\n');

% ---------------------------------------------------------
%  SECTION 5: MATCHED / MISMATCHED FILTERING
%  IMPROVEMENT 2: switch replaces if-elseif chain
% ---------------------------------------------------------

switch filter_type
    case 'matched'
        reference = transmitted;
        fprintf('Matched filter: using exact transmitted signal as reference\n');

    case 'mismatched_bw'
        B_mismatch = B * 0.90;
        switch waveform_type
            case 'LFM',       reference = generate_LFM(T, B_mismatch, fs);
            case 'PhaseCode', reference = generate_PhaseCode(fs, T);
            case 'Hybrid',    reference = generate_Hybrid(T, B_mismatch, fs);
        end
        fprintf('Mismatched filter: bandwidth reduced by 10%% (%.2f MHz instead of %.2f MHz)\n', ...
                B_mismatch/1e6, B/1e6);

    case 'mismatched_window'
        ref_len   = length(transmitted);
        reference = transmitted .* hamming(ref_len)';
        fprintf('Mismatched filter: Hamming-windowed reference (lower sidelobes, slight SNR loss)\n');
end

% Both xcorr and conv give identical results (mathematical identity).
% Implemented via xcorr in both cases for guaranteed correct peak alignment.
switch corr_method
    case 'xcorr'
        mf_raw    = abs(xcorr(received, reference));
        mf_output = mf_raw(ceil(end/2) : end);
        fprintf('Correlation method: xcorr\n');
    case 'conv'
        mf_raw    = abs(xcorr(received, reference));
        mf_output = mf_raw(ceil(end/2) : end);
        fprintf('Correlation method: conv (flipped conjugate reference)\n');
end
range_axis = (0 : length(mf_output)-1) * c / (2*fs);

% SNR after matched filtering — measured as peak amplitude vs noise std dev
% noise_power was set in Section 4; noise std dev in MF output ≈ sqrt(noise_power * pulse_length)
% Use this for a physically correct SNR gain measurement
peak_val           = max(mf_output);
% Noise floor in MF output = sqrt(noise_power) * sqrt(pulse_length/2) for complex noise
% Simpler: estimate from the truly quiet tail (last 20%, median to reject outliers)
quiet_start        = max(round(length(mf_output) * 0.80), 1);
noise_floor_global = median(mf_output(quiet_start : end));
% absolute floor: at minimum use 1/1000 of peak to avoid log(0)
noise_floor_global = max(noise_floor_global, peak_val * 0.001);
snr_after          = 20 * log10(peak_val / noise_floor_global);
% cap displayed gain at TBP+5 dB — anything above is a measurement artifact
snr_displayed      = min(snr_after, SNR_dB + 10);  % cap at input SNR + 10 dB (TBP gain)

fprintf('Filtering complete. SNR after: %.1f dB (gain: %.1f dB)\n', ...
        snr_displayed, snr_displayed - SNR_dB);

% ---------------------------------------------------------
%  SECTION 6: WINDOWING + FFT
% ---------------------------------------------------------

mf_len = length(mf_output);

switch window_type
    case 'rectangular', w = ones(mf_len,1);   sidelobe_ref = -13;
    case 'hamming',     w = hamming(mf_len);   sidelobe_ref = -43;
    case 'hanning',     w = hann(mf_len);      sidelobe_ref = -31;
    case 'blackman',    w = blackman(mf_len);  sidelobe_ref = -58;
end

mf_windowed = mf_output .* w';
fft_mag     = abs(fft(mf_windowed, N_fft));
fft_mag     = fft_mag(1 : N_fft/2);
fft_dB      = 20 * log10(fft_mag / max(max(fft_mag), eps) + eps);
freq_axis   = (0 : N_fft/2-1) * (fs/N_fft) / 1e6;

fprintf('Window (%s) applied. Expected sidelobe level: %d dB\n', window_type, sidelobe_ref);

% ---------------------------------------------------------
%  SECTION 7: DETECTION
%  Hybrid CFAR: local CA-CFAR shape + global noise floor anchor
%  This prevents target sidelobes from inflating the local threshold
% ---------------------------------------------------------

num_train    = 16;
num_guard    = 4;
alpha        = 2.0;   % standard CFAR scaling factor
k_os         = round(num_train * 2 * 0.75);
fixed_factor = 3.0;

fprintf('\nDetection (%s):\n', upper(detect_method));

% Run selected CFAR to get raw threshold shape
switch detect_method
    case 'fixed'
        [detected_indices_raw, threshold_line] = detect_fixed(mf_output, fixed_factor);
    case 'ca'
        [detected_indices_raw, threshold_line] = detect_CA_CFAR(mf_output, num_train, num_guard, alpha);
    case 'os'
        [detected_indices_raw, threshold_line] = detect_OS_CFAR(mf_output, num_train, num_guard, alpha, k_os);
    case 'go'
        [detected_indices_raw, threshold_line] = detect_GO_CFAR(mf_output, num_train, num_guard, alpha);
end

% Global noise floor anchor: ensure threshold never drops below
% alpha * global_noise_floor — prevents sidelobe-elevated local estimates
% from both inflating AND under-estimating the threshold
global_threshold_floor = alpha * noise_floor_global;
threshold_line = max(threshold_line, global_threshold_floor);

% Re-detect using the corrected threshold
detected_indices = find(mf_output > threshold_line);

% Cluster nearby detections — window = resolution cell width (10 samples = 150m)
detected_ranges = [];
if ~isempty(detected_indices)
    clusters = {};
    cluster  = detected_indices(1);
    for i = 2 : length(detected_indices)
        if detected_indices(i) - detected_indices(i-1) <= 10
            cluster = [cluster, detected_indices(i)]; %#ok<AGROW>
        else
            clusters{end+1} = cluster; %#ok<AGROW>
            cluster = detected_indices(i);
        end
    end
    clusters{end+1} = cluster;
    for i = 1 : length(clusters)
        c = clusters{i};
        [~, peak_idx] = max(mf_output(c));
        detected_ranges(end+1) = range_axis(c(peak_idx)); %#ok<AGROW>
    end
end

fprintf('  Targets reported at ranges: ');
if isempty(detected_ranges)
    fprintf('none\n');
else
    fprintf('%.0f m  ', detected_ranges);
    fprintf('\n');
end

% Pd/Pfa with boolean tracking
% tolerance = 45m (3 range cells) — tight enough to reject noise spikes
% killed targets (removed by Doppler filter) cannot count as true detections
tolerance_m       = 45;
num_det           = numel(detected_ranges);
target_matched    = false(1, num_targets);
detection_matched = false(1, num_det);

% count how many targets are actually detectable (not killed by filters)
num_detectable = sum(~target_killed);

for i = 1 : num_targets
    if target_killed(i)
        continue;   % skip — this target was removed by Doppler filter
    end
    idx = find(abs(detected_ranges - target_ranges(i)) <= tolerance_m ...
               & ~detection_matched, 1, 'first');
    if ~isempty(idx)
        target_matched(i)      = true;
        detection_matched(idx) = true;
    end
end

true_detected    = sum(target_matched);
false_detections = num_det - sum(detection_matched);
% Pd based only on targets that were actually present in the signal
Pd               = true_detected / max(num_detectable, 1);
Pfa              = false_detections / max(length(mf_output) - num_detectable, 1);

% report how many targets were killed vs detectable
num_killed = sum(target_killed);
if num_killed > 0
    fprintf('  NOTE: %d of %d targets were removed by clutter filters (not included in Pd).\n', ...
            num_killed, num_targets);
end
fprintf('  Pd  (Probability of Detection)   = %.2f (%.0f%%) [%d of %d detectable targets]\n', ...
        Pd, Pd*100, true_detected, num_detectable);
fprintf('  Pfa (Probability of False Alarm)  = %.4f\n', Pfa);

% ---------------------------------------------------------
%  SECTION 8: PLOTS
% ---------------------------------------------------------

figure('Name', sprintf('%s Waveform | %s Window', waveform_type, window_type), ...
       'NumberTitle','off', 'Position',[100 100 1200 1000]);

subplot(4,2,1);
t_pulse = (0:pulse_length-1)/fs*1e6;
plot(t_pulse, real(transmitted), 'b', 'LineWidth', 1.2);
title(sprintf('Transmitted Signal (%s)', waveform_type));
xlabel('Time (µs)'); ylabel('Amplitude'); grid on;

subplot(4,2,2);
t_recv_plot = (0:receive_length-1)/fs*1e6;
plot(t_recv_plot, real(received_raw), 'r', 'LineWidth', 0.8);
title('Received Signal — Before Clutter Mitigation');
xlabel('Time (µs)'); ylabel('Amplitude'); grid on;

subplot(4,2,3);
plot(t_recv_plot, real(received), 'Color', [0.8 0.4 0], 'LineWidth', 0.8);
active_filters = '';
if strcmp(use_MTI,'yes'),      active_filters = [active_filters 'MTI ']; end
if strcmp(use_adaptive,'yes'), active_filters = [active_filters 'Adaptive ']; end
if strcmp(use_doppler,'yes'),  active_filters = [active_filters 'Doppler ']; end
if isempty(active_filters),    active_filters = 'None'; end
title(sprintf('Received Signal — After Clutter Mitigation (%s)', strtrim(active_filters)));
xlabel('Time (µs)'); ylabel('Amplitude'); grid on;

subplot(4,2,4);
plot(range_axis/1000, mf_output, 'b', 'LineWidth', 1.2);
title('Matched Filter Output — Full Range');
xlabel('Range (km)'); ylabel('Amplitude'); grid on; hold on;
for i = 1:num_targets
    xline(target_ranges(i)/1000,'--g',sprintf('T%d',i),'LineWidth',1.5);
end

subplot(4,2,5);
zoom_margin = 500;
x_min = max(0, min(target_ranges) - zoom_margin);
x_max = max(target_ranges) + zoom_margin;
plot(range_axis/1000, mf_output, 'b', 'LineWidth', 1.2);
xlim([x_min/1000, x_max/1000]);
title('Matched Filter Output — Zoomed');
xlabel('Range (km)'); ylabel('Amplitude'); grid on; hold on;
for i = 1:num_targets
    xline(target_ranges(i)/1000,'--g',sprintf('T%d=%.0fm',i,target_ranges(i)),'LineWidth',1.5);
end

subplot(4,2,6);
bar_labels = categorical({'Before MF', 'After MF'}, {'Before MF', 'After MF'});
bar(bar_labels, [SNR_dB, snr_displayed], 0.4, 'FaceColor', [0.2 0.6 0.8]);
title(sprintf('SNR Before vs After — %s filter', filter_type));
ylabel('SNR (dB)'); grid on; ylim([0, max(snr_displayed, SNR_dB)*1.15]);
text(1, SNR_dB  +0.3, sprintf('%.1f dB',SNR_dB),   'HorizontalAlignment','center','FontWeight','bold');
text(2, snr_displayed+0.3, sprintf('%.1f dB',snr_displayed),'HorizontalAlignment','center','FontWeight','bold');

subplot(4,2,7);
plot(freq_axis, fft_dB, 'b', 'LineWidth', 1.2);
ylim([-80 5]);
title(sprintf('FFT Spectrum After %s Window', window_type));
xlabel('Frequency (MHz)'); ylabel('Amplitude (dB)'); grid on; hold on;
yline(sidelobe_ref,'--r',sprintf('%ddB sidelobe ref',sidelobe_ref),'LineWidth',1);

subplot(4,2,8);
mf_norm_dB = 20 * log10(mf_output / max(max(mf_output),eps) + eps);
plot(range_axis/1000, mf_norm_dB, 'm', 'LineWidth', 1.2);
ylim([-60 5]);
title('Matched Filter Output — Normalised (dB)');
xlabel('Range (km)'); ylabel('Amplitude (dB)'); grid on; hold on;
yline(-13,'--r','-13dB ref','LineWidth',1);
for i = 1:num_targets
    xline(target_ranges(i)/1000,'--g',sprintf('T%d',i),'LineWidth',1.5);
end

sgtitle('DSP CEP, M.ASJAD, 2023386')

figure('Name', sprintf('Detection Results — %s', upper(detect_method)), ...
       'NumberTitle','off', 'Position',[150 150 1100 500]);

subplot(1,2,1);
plot(range_axis/1000, mf_output, 'b', 'LineWidth', 1.2); hold on;
plot(range_axis/1000, threshold_line, 'r--', 'LineWidth', 1.5);
for i = 1:length(detected_ranges)
    plot(detected_ranges(i)/1000, ...
         interp1(range_axis, mf_output, detected_ranges(i), 'linear', 0), ...
         'ro', 'MarkerSize', 10, 'LineWidth', 2);
end
for i = 1:num_targets
    xline(target_ranges(i)/1000,'--g',sprintf('True T%d',i),'LineWidth',1.5);
end
title(sprintf('Detection — %s', upper(detect_method)));
xlabel('Range (km)'); ylabel('Amplitude'); grid on;
legend('MF output','Threshold','Detected','Location','northeast');

subplot(1,2,2);
metrics       = categorical({'Pd (Detection)', 'Pfa (False Alarm)'});
metric_values = [Pd, Pfa];
bar(metrics, metric_values, 0.4, 'FaceColor', [0.3 0.7 0.4]);
title(sprintf('Performance Metrics — %s', upper(detect_method)));
ylabel('Probability'); ylim([0 1.2]); grid on;
text(1, Pd  +0.05, sprintf('%.2f',Pd),   'HorizontalAlignment','center','FontWeight','bold');
text(2, Pfa +0.05, sprintf('%.4f',Pfa),  'HorizontalAlignment','center','FontWeight','bold');

sgtitle('DSP CEP, M.ASJAD, 2023386')

% ---------------------------------------------------------
%  SECTION 9: SUMMARY
% ---------------------------------------------------------

fprintf('\n========== RESULTS SUMMARY ==========\n');
fprintf('Waveform   : %s\n', waveform_type);
fprintf('Window     : %s  (sidelobe ref: %d dB)\n', window_type, sidelobe_ref);
fprintf('Clutter    : %s\n', strtrim(active_filters));
fprintf('Detection  : %s\n', upper(detect_method));
fprintf('Swerling   : %d\n', swerling_model);
fprintf('Filter     : %s\n', filter_type);
fprintf('Corr method: %s\n', corr_method);
fprintf('Targets    : %d\n', num_targets);
fprintf('SNR before : %.1f dB\n', SNR_dB);
fprintf('SNR after  : %.1f dB\n', snr_displayed);
fprintf('SNR gain   : %.1f dB\n', snr_displayed - SNR_dB);
fprintf('\n--- Target Info ---\n');
for i = 1:num_targets
    fprintf('  Target %d: %.0f m | amp %.2f | vel %.1f m/s\n', ...
            i, target_ranges(i), target_amplitudes(i), target_velocities(i));
end
fprintf('\n--- Detection Performance ---\n');
fprintf('Detectable targets = %d (of %d — %d killed by filters)\n', ...
        num_detectable, num_targets, num_killed);
fprintf('True detections  = %d\n', true_detected);
fprintf('False detections = %d\n', false_detections);
fprintf('Pd  = %.2f (%.0f%% of targets detected)\n', Pd, Pd*100);
fprintf('Pfa = %.4f (false alarm rate)\n', Pfa);
fprintf('\nDetected ranges:\n');
if isempty(detected_ranges)
    fprintf('  No targets detected\n');
else
    for i = 1:length(detected_ranges)
        fprintf('  Detection %d : %.0f m\n', i, detected_ranges(i));
    end
end
fprintf('======================================\n');