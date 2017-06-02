function figure_wire(surf, Qe, plain)
   
    %figure_wire(surf, Ecolor, Fcolor);
    %
    % The function colors mesh faces and edges 
    %
    % surf  : surface mesh
    % Ecolor: edge color 
    % Fcolor: face color
    %
    % For instance, we can have Ecolor=[0.8 0.8 0.8]
    %
    % (C) 2008 Moo K. Chung 
    %  mkchung@wisc.edu
    %  Department of Biostatisics and Medical Informatics
    %  University of Wisconsin, Madison
    % 
    % 2008 created
    % 2013 Sept. 5 nargin added

    background = 'white';
    whitebg(gcf,background);

    p = patch(surf);
    if ~plain
        set(p,'CData',Qe,'FaceColor','interp','EdgeColor','none');
    else
        set(p,'CData',Qe,'FaceColor','w','EdgeColor','b');
    end
    caxis([0,30]);

    daspect([1 1 1]); 
    axis off;

    set(gcf, 'Color', background, 'InvertHardcopy', 'off');

    camlight;
    lighting gouraud;
    material shiny; 
end