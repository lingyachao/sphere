clear; close all;

load('./computed_brain_grid/N40962.mat');
avg_D = 0.3;

% set the output directory
folder_name = ['brain_N' num2str(N) '_' datestr(now, 'mmddHHMM')];
OUTPUT_DIR = ['./data/' folder_name '/'];
mkdir(OUTPUT_DIR);

% initialize parameters, initial state and map
k = 0;
K = 200;
T0 = 1;
visualize = true;
map = make_map(laplacian);
last = make_IC(N);

% set plotting window
f = figure;
set(f, 'Position', [200 0 600 800]);
load('autism.surface.mat', 'tri');
surf.vertices = locs;
surf.faces = tri;

for k = 1:K
    
    if k < 0 / T0
        source_drive = mean(last.dVe(:));
    elseif k > 150 / T0
        source_drive = mean(last.dVe(:));
    else
        source_drive = 3;
    end

    fprintf(['Running simulation , ' num2str(k) ' ... \n']);
    tic; % start timer
    [time,last,Qe_micro,Qe_macro] = seizing_cortical_field(...
        source_drive, map, T0, last, ...
        laplacian, avg_D, [], []);
    fprintf(['Run time ' num2str(toc) '\n']);

    if visualize
        figure_wire(surf, last.Qe);
        drawnow;
    end
    
    Qe = last.Qe;
    save([OUTPUT_DIR 'seizing_cortical_field_k_' num2str(k) '.mat'], ...
        'time', 'Qe', 'Qe_micro', 'Qe_macro');
end 
