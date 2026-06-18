function output = apply_doppler(received, fs)
% Doppler-Based Discrimination Filter
% Zeros out zero-Doppler region in frequency domain where clutter sits
%
% FIX: notch is now defined in Hz not as a fraction of fs
% Old method: 2% of fs = 200,000 Hz notch — this killed ALL moving targets
% New method: fixed 2000 Hz notch — kills stationary clutter only
% Targets moving faster than ~1.5 m/s survive the notch

    N            = length(received);
    fft_received = fft(received, N);

    % notch width in Hz — kills clutter within ±1000 Hz of zero Doppler
    % corresponds to targets slower than ~0.15 m/s (effectively stationary)
    notch_hz   = 2000;                       % total notch width in Hz
    notch_bins = round(notch_hz * N / fs);   % convert Hz to number of FFT bins
    notch_bins = max(notch_bins, 1);         % at least 1 bin

    fft_received(1 : notch_bins)       = 0;
    fft_received(N-notch_bins+1 : end) = 0;

    output = ifft(fft_received, N);

    fprintf('  Doppler filter applied (notch: %.0f Hz = %d bins)\n', notch_hz, notch_bins);
end
