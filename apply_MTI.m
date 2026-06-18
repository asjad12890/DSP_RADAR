function output = apply_MTI(received)
% MTI (Moving Target Indicator) Filter
% Subtracts consecutive samples to cancel stationary clutter
% Uses a 3-pulse canceller for stronger clutter suppression

    N      = length(received);
    output = zeros(1, N);

    for n = 3 : N
        output(n) = received(n) - 2*received(n-1) + received(n-2);
    end

    fprintf('  MTI filter applied (3-pulse canceller)\n');
end
