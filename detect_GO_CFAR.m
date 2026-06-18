function [detected_indices, threshold_line] = detect_GO_CFAR(mf_output, num_train, num_guard, alpha)
% GO-CFAR: takes greater of left and right training cell averages

    N                = length(mf_output);
    threshold_line   = zeros(1, N);
    detected_indices = [];
    half_window      = num_train + num_guard;

    for cut = half_window+1 : N-half_window
        left_train  = mf_output(cut - half_window : cut - num_guard - 1);
        right_train = mf_output(cut + num_guard + 1 : cut + half_window);
        noise_estimate = max(mean(left_train), mean(right_train));
        threshold_line(cut) = alpha * noise_estimate;
        if mf_output(cut) > threshold_line(cut)
            detected_indices = [detected_indices, cut];
        end
    end

    % suppressed per-call print
end