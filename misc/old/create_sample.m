function [data_macro,data_micro,fs] = create_sample(DATA_DIR)

    load([DATA_DIR 'seizing_cortical_field_k_1.mat']);

    % starting k for each period
    k0 = 10:10:140;
    % number of files for each period
    K = 10;
    % number of time points in each file
    T = size(Qe_macro, 1);
    % number of nodes
    N_macro = size(Qe_macro, 2);
    N_micro = size(Qe_micro, 2);

    dec = 10;
    data_macro = zeros(K*T/10, N_macro, length(k0));
    data_micro = zeros(K*T/10, N_micro, length(k0));

    macro = zeros(K*T, N_macro);
    micro = zeros(K*T, N_micro);

    for i = 1:length(k0)

        for k = 1:K
            fprintf(['Read in ' num2str(k0(i) + k) '\n'])
            load([DATA_DIR 'seizing_cortical_field_k_'  num2str(k0(i) + k) '.mat'])
            macro(1+(k-1)*T : k*T,:) = Qe_macro;
            micro(1+(k-1)*T : k*T,:) = Qe_micro;
        end

        for k = 1:N_macro
            data_macro(:,k,i) = decimate(macro(:,k), dec);
        end

        for k = 1:N_micro
            data_micro(:,k,i) = decimate(micro(:,k), dec);
        end
        
    end

    % sampling interval
    dt = (time(10) - time(9)) * dec;
    % sampling frequency
    fs = 1/dt;

    % save the results
    save([DATA_DIR 'downsample_data.mat'], 'data_macro', 'data_micro', 'fs')
end