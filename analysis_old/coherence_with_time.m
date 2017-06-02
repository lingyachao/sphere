function [avg_coh_macro, avg_coh_micro] = coherence_with_time(data_macro, data_micro, fs)

    % set up for computing coherence
    CHRONUX_PATH = './chronux_2_12'; 
    addpath(genpath(CHRONUX_PATH));

    BAND = [1 13];                  % select a frequency range to analyze
    TW = 20;                        % time-bandwidth product
    ntapers = 2*TW-1;               % choose the # of tapers.
    params.tapers = [TW, ntapers];  % ... time-bandwidth product and tapers.
    params.Fs = fs;                 % ... sampling rate
    params.pad = -1;                % ... no zero padding.
    params.fpass = BAND;            % ... freq range to pass
    params.err = [1 0.05];          % ... theoretical error bars, p=0.05.

    avg_coh_macro = zeros(size(data_macro,3), 1);
    avg_coh_micro = zeros(size(data_micro,3), 1);

    for i = 1:size(data_macro,3)
        [coh,~,~,~] = compute_coherence(data_macro(:,:,i), params);
        avg_coh_macro(i) = mean(mean(mean(coh)));
        [coh,~,~,~] = compute_coherence(data_micro(:,:,i), params);
        avg_coh_micro(i) = mean(mean(mean(coh)));
    end
    
end