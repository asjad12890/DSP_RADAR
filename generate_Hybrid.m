function signal = generate_Hybrid(T, B, fs)
% Generates a Hybrid waveform: LFM chirp modulated by Barker-13 phase code
% Each chip is a mini LFM that sweeps bandwidth B over its chip duration T_chip
% This gives LFM-like range resolution (c/2B) with PhaseCode-like sidelobes

    barker13         = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1];
    num_chips        = length(barker13);
    total_samples    = round(T * fs);
    samples_per_chip = floor(total_samples / num_chips);
    T_chip           = samples_per_chip / fs;

    % Each chip sweeps full bandwidth B over its duration T_chip
    % chirp rate per chip = B / T_chip  (much faster than B/T)
    k_chip = B / T_chip;

    signal = zeros(1, total_samples);

    for i = 1 : num_chips
        start_idx  = (i-1) * samples_per_chip + 1;
        end_idx    = min(i * samples_per_chip, total_samples);
        n_samples  = end_idx - start_idx + 1;
        t_chip     = (0 : n_samples-1) / fs;
        chip       = barker13(i) * exp(1j * pi * k_chip .* t_chip.^2);
        signal(start_idx : end_idx) = chip;
    end

    signal = signal / max(abs(signal));
end