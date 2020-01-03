load('Y:\Projects\Seizures\Modeling_Seizure_Dynamics_(Monica)\sphere\data\sphere_N10242_R10_12211753_Nie670_source8for10s\data_sample.mat');
Qe_tracing = Qe_rand(:,2);

fine_time = (1:150000)/500;
[pks,locs,w,p] = findpeaks(Qe_tracing, fine_time, 'MinPeakHeight', 5);

ictal_locs = locs(locs > 150);
intervals = ictal_locs(2:end) - ictal_locs(1:end-1);

figure;
subplot(1,2,1);
plot(intervals);

ictal_widths = w(locs > 150);
subplot(1,2,2);
plot(ictal_widths);
