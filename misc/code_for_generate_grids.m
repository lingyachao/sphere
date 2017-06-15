%% generate sphere grid

[N,locs,laplacian,avg_D] = make_sphere(10000, R);
save(['./computed_sphere_grid/N' num2str(N) '_R' num2str(R) '.mat'], ...
    'N', 'locs', 'laplacian', 'avg_D', 'R');

%% generate sphere grid with hole

[N,locs,laplacian,avg_D] = make_sphere_hole(10000, R);
save(['./computed_sphere_grid/N' num2str(N) '_R' num2str(R) '_hole.mat'], ...
    'N', 'locs', 'laplacian', 'avg_D', 'R', 'micro_idx', 'macro_idx');
load('./computed_sphere_grid/N9747_R10_hole.mat');
micro_idx = [5121,5311,5312,5120,5122,4931,4932]';
macro_idx = [5121,5803,5804,5117,5125,4439,4440,5493,5494,5119,5123,4749,4750]';

%% generate brain grid

[N,locs,laplacian,avg_D] = make_brain_grid;
save(['./computed_brain_grid/N' num2str(N) '.mat'], ...
    'N', 'locs', 'laplacian', 'avg_D', ...
    'micro_idx', 'macro_idx', 'focus_idx', 'normal_idx');

%% generate small sphere
R = 10;
[N,locs,laplacian,avg_D] = make_sphere(100, R);
micro_idx = [];
macro_idx = [];
focus_idx = 1;
normal_idx = 40;
save(['./computed_sphere_grid/N' num2str(N) '_R' num2str(R) '.mat'], ...
    'N', 'locs', 'laplacian', 'avg_D', 'R', ...
    'micro_idx', 'macro_idx', 'focus_idx', 'normal_idx');
