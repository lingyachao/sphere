surf.vertices=squeeze(controlouter(1,:,:))';
surf.faces=tri;
figure; figure_wire(surf, 'k', 'w');

hold on;
for i = 1 : 10
    scatter3(controlouter(1,1,i), controlouter(1,2,i), controlouter(1,3,i), 30, 'r', 'filled');
end

a = NaN(40962, 1);

for i = 1 : 40962
    a(i) = sum(sum(tri == i));
end