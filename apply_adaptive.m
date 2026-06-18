function output = apply_adaptive(received)
% Adaptive Filter using LMS (Least Mean Squares) algorithm
% Learns what the clutter looks like then subtracts it
% Works well against non-stationary clutter that MTI cannot handle
%
% FIX: mu reduced from 0.01 to 0.001 to prevent target cancellation
% Smaller mu = slower adaptation = less risk of cancelling target echo

    N          = length(received);
    filter_len = 8;       % how many past samples the filter looks at
    mu         = 0.001;   % FIXED: was 0.01 — too aggressive, was eating target signal
                          % 0.001 adapts slowly but safely preserves target echoes

    weights    = zeros(1, filter_len);
    output     = zeros(1, N);

    for n = filter_len : N
        x                = received(n : -1 : n-filter_len+1);
        clutter_estimate = weights * x.';
        error_signal     = received(n) - clutter_estimate;
        weights          = weights + mu * error_signal * conj(x);
        output(n)        = error_signal;
    end

    fprintf('  Adaptive filter applied (LMS, filter_len=%d, mu=%.4f)\n', filter_len, mu);
end
