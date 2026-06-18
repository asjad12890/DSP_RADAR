function signal = generate_PhaseCode(fs, T)
% Generates a phase-coded waveform using Barker-13 code

    barker13         = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1];
    num_chips        = length(barker13);
    total_samples    = round(T * fs);
    samples_per_chip = floor(total_samples / num_chips);

    signal = zeros(1, total_samples);

    for i = 1 : num_chips
        start_idx = (i-1) * samples_per_chip + 1;
        end_idx   = min(i * samples_per_chip, total_samples);
        signal(start_idx : end_idx) = barker13(i);
    end

    signal = signal / max(abs(signal));
end
