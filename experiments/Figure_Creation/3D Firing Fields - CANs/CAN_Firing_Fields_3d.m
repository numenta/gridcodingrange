

data = load("./data/data1.mat");

f = data.arr;
figure; hold on;

thresh = 0.012


b = data.b


x = -1:2./(b-1.):1.;

y = x;
z = x;

[x,y,z] = meshgrid(x,y,z);

size(x)
size(f)

% % 
set(gca,'linewidth',4)
set(gca, 'XTickLabel', [])
set(gca, 'YTickLabel', [])
set(gca, 'ZTickLabel', [])
% 
xticks([-1 0 1])
yticks([-1 0 1])
zticks([-1 0 1])

subplot(1,1,1);


isosurface(x, y, z, f, thresh);
isocaps(x, y, z, f, thresh, 'above');



azim = -25;
elev = 30;
lightangle(azim,elev);
daspect([1,1,1]);
camlight
lighting gouraud
grid on




view(22.5,15)

% axis vis3d
colormap jet
% colormap copper



v1 = [1 1 -1;1 -1 -1; -1 -1 -1; -1 1 -1];
f1 = [1 2 3 4];
patch('Faces',f1,'Vertices',v1,'FaceColor','black','FaceAlpha',.15, 'LineStyle','none');

v1 = [1  1 -1;
      1  1  1; 
      -1 1 1; 
      -1 1 -1];
f1 = [1 2 3 4];
patch('Faces',f1,'Vertices',v1,'FaceColor','black','FaceAlpha',.05, 'LineStyle','none');


v1 = [-1  -1 -1;
      -1   1  -1; 
      -1  1  1; 
      -1 -1 1];
f1 = [1 2 3 4];
patch('Faces',f1,'Vertices',v1,'FaceColor','black','FaceAlpha',.05, 'LineStyle','none');



v1 = [1 1 1;1 -1 1; -1 -1 1; -1 1 1];
f1 = [1 2 3 4];
patch('Faces',f1,'Vertices',v1,'FaceColor','black','FaceAlpha',0., 'LineWidth',4.);

v1 = [1 -1 1;1 -1 -1;-1 -1  -1;-1 -1  1];
f1 = [1 2 3 4];
patch('Faces',f1,'Vertices',v1,'FaceColor','None','FaceAlpha',1., 'LineWidth',4.);

v1 = [1 1 1 ; 
      1 1 -1];
f1 = [1 2];
patch('Faces',f1,'Vertices',v1,'FaceColor','None','FaceAlpha',1., 'LineWidth',4.);

v1 = [1 -1  -1 ; 
      1  1 -1];
f1 = [1 2];
patch('Faces',f1,'Vertices',v1,'FaceColor','None','FaceAlpha',1., 'LineWidth',4.);




M = max(f,[],3);

planeimg = mat2gray(reshape(M,100,100));
size(planeimg)


% plot the image plane using surf.
surf([-1 1], [-1 1], repmat(-1, [2 2]), planeimg,'facecolor','texture')



saveas(gcf, './Figures/3d_field.png');


