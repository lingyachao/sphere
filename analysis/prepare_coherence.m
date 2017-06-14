% required toolbox
% CHRONUX_PATH = './chronux_2_12';
% addpath(genpath(CHRONUX_PATH));

% parameters for computing coherence
BAND = [1 13];                   % select a frequency range to analyze
TW = 20;                         % time-bandwidth product
ntapers = 2*TW-1;                % choose the # of tapers.
params.tapers = [TW, ntapers];   % ... time-bandwidth product and tapers.
params.Fs = double(500);         % ... sampling rate
params.pad = -1;                 % ... no zero padding.
params.fpass = BAND;             % ... freq range to pass
params.err = [1 0.05];           % ... theoretical error bars, p=0.05.

% allocate space to store coherence data
[~,~,freq,~] = compute_coherence(Qe_macro(1:per_P,1:2), params);
N_freq = length(freq); 

% subset macro/micronodes to reduce run time
% Qe_macro = Qe_macro(:,[1:7,8:2:18]);
% Qe_micro = Qe_micro(:,[1:7,8:2:18]);

N_macro = size(Qe_macro, 2);
macro_t_coh = NaN(N_macro, N_macro, N_freq, P);
macro_t_coh_conf = NaN(N_macro, N_macro, N_freq, P);
macro_t_phi = NaN(N_macro, N_macro, N_freq, P);

N_micro = size(Qe_micro, 2);
micro_t_coh = NaN(N_micro, N_micro, N_freq, P);
micro_t_coh_conf = NaN(N_micro, N_micro, N_freq, P);
micro_t_phi = NaN(N_micro, N_micro, N_freq, P);

for p = 1:P
    fprintf(['Computing coherence for period ' num2str(p) '\n']);
    % range = ((p-1) * per_P / 2) + 1 : ((p+1) * per_P / 2);
    range = ((p-1) * per_P) + 1 : (p*per_P);
    [coh,phi,~,coh_conf] = compute_coherence(Qe_macro(range,:), params);
    macro_t_coh(:,:,:,p) = coh;
    macro_t_coh_conf(:,:,:,p) = coh_conf;
    macro_t_phi(:,:,:,p) = phi;
    
    [coh,phi,~,coh_conf] = compute_coherence(Qe_micro(range,:), params);
    micro_t_coh(:,:,:,p) = coh;
    micro_t_coh_conf(:,:,:,p) = coh_conf;
    micro_t_phi(:,:,:,p) = phi;
end

save(COHERENCE_FILE, 'macro_t_coh', 'macro_t_coh_conf', 'macro_t_phi', 'freq', ...
                     'micro_t_coh', 'micro_t_coh_conf', 'micro_t_phi');