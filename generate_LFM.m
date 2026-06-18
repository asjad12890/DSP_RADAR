function signal = generate_LFM(T, B, fs)
% Generates a Linear Frequency Modulated (LFM) chirp signal
    t      = 0 : 1/fs : T - 1/fs;
    k      = B / T;
    signal = exp(1j * pi * k .* t.^2);
end
