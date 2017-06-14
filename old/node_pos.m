r1 = 1;
r2 = 2;

x_center = 0;
y_center = 0;
thetas = pi/2 : pi/3 : pi/2 + 5*pi/3;
x_level1 = r1 * cos(thetas);
y_level1 = r1 * sin(thetas);

x_level2_half = r2 * cos(thetas);
y_level2_half = r2 * sin(thetas);
x_level2 = nan(1,12);
y_level2 = nan(1,12);
x_level2(1:2:11) = x_level2_half;
y_level2(1:2:11) = y_level2_half;
x_level2(2:2:12) = (x_level2_half + x_level2_half([2:6,1])) / 2;
y_level2(2:2:12) = (y_level2_half + y_level2_half([2:6,1])) / 2;

x = [x_center, x_level1, x_level2];
y = [y_center, y_level1, y_level2];
position = [x;y]';

save('./computed_sphere_grid/node_pos_19.mat', 'position');