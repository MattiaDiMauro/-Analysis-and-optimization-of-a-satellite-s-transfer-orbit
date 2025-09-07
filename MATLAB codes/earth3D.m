function earth3D

image_file = 'http://upload.wikimedia.org/wikipedia/commons/thumb/c/cd/Land_ocean_ice_2048.jpg/1024px-Land_ocean_ice_2048.jpg';


[X,Y,Z] = sphere;
X = X * earthRadius * 1e-3;
Y = Y * earthRadius * 1e-3;
Z = Z * earthRadius * 1e-3;

globe = surf(X, Y, -Z);
hold on
grid on
view(135,15);
alpha = 1;

cdata = imread(image_file);
set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');


% axis of cartesian reference system:
quiver3(zeros(3,1),zeros(3,1),zeros(3,1),[2*earthRadius*1e-3;0;0],...
    [0;2*earthRadius*1e-3;0],[0;0;2*earthRadius*1e-3],'LineWidth',1.15,...
    'Color','blue');

% name of axis:
text(2*earthRadius*1e-3+120,100,-100,texlabel('gamma'),'Color','k','FontSize',15);
text(100,2*earthRadius*1e-3+120,-100,texlabel('Y'),'Color','k','FontSize',10);
text(100,100,2*earthRadius*1e-3+140,texlabel('Z'),'Color','k','FontSize',10);

axis equal
end