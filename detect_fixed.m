function [detected_indices, threshold_line] = detect_fixed(mf_output, threshold_factor)
% Fixed Threshold Detection — baseline comparison method

    % FIX: use last 10% of signal as noise estimate (after all echoes)
    % this avoids the clutter-contaminated early samples
    noise_region = mf_output(round(end*0.85) : end);
    noise_floor  = mean(noise_region);

    threshold      = threshold_factor * noise_floor;
    threshold_line = threshold * ones(1, length(mf_output));
    detected_indices = find(mf_output > threshold);

    % suppressed per-call prints
end