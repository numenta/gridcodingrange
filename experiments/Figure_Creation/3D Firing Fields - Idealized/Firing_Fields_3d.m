
x = -1:0.02:1.;
y = x;
z = x;
[x,y,z] = meshgrid(x,y,z);
sz = size(x);

x_ = reshape(x,1,[]);
y_ = reshape(y,1,[]);
z_ = reshape(z,1,[]);
X = vertcat(x_,y_,z_);


m = 2
f = get_bubbles(m, X);
f = reshape(f, sz);

R = 0.9;
min_f=min(f, [],'all');
max_f=max(f, [],'all');
thresh=(1 - R)*min_f + R*max_f;



set(gca,'linewidth',4)
set(gca, 'XTickLabel', [])
set(gca, 'YTickLabel', [])
set(gca, 'ZTickLabel', [])
xticks([-1 0 1])
yticks([-1 0 1])
zticks([-1 0 1])


subplot(1,1,1);
% -------------
isosurface(x, y, z, f, thresh);
isocaps(x, y, z, f, thresh, 'above');

xlim([-1 1]);
ylim([-1 1]);
zlim([-1 1]);

azim = -25;
elev = 30;
lightangle(azim,elev);
daspect([1,1,1]);
camlight
lighting gouraud
grid on
view(22.5,15)
colormap jet

v1 = [1 1 -1;1 -1 -1; -1 -1 -1; -1 1 -1];
f1 = [1 2 3 4];
patch('Faces',f1,'Vertices',v1,'FaceColor','black','FaceAlpha',.15, 'LineStyle','none');
v1 = [1  1 -1;
      1  1  1; 
      -1 1  1; 
      -1 1 -1];
f1 = [1 2 3 4];
patch('Faces',f1,'Vertices',v1,'FaceColor','black','FaceAlpha',.05, 'LineStyle','none');

v1 = [-1  -1 -1;
      -1   1 -1; 
      -1   1  1; 
      -1  -1  1];
f1 = [1 2 3 4];
patch('Faces',f1,'Vertices',v1,'FaceColor','black','FaceAlpha',.05, 'LineStyle','none');

v1 = [1 1 1;1 -1 1; -1 -1 1; -1 1 1];
f1 = [1 2 3 4];
patch('Faces',f1,'Vertices',v1,'FaceColor','black','FaceAlpha',0., 'LineWidth',4.);

v1 = [1 -1 1;1 -1 -1;-1 -1  -1;-1 -1  1];
f1 = [1 2 3 4];
patch('Faces',f1,'Vertices',v1,'FaceColor','None','FaceAlpha',1., 'LineWidth',4.);

v1 = [1 1  1; 
      1 1 -1];
f1 = [1 2];
patch('Faces',f1,'Vertices',v1,'FaceColor','None','FaceAlpha',1., 'LineWidth',4.);

v1 = [1 -1 -1; 
      1  1 -1];
f1 = [1 2];
patch('Faces',f1,'Vertices',v1,'FaceColor','None','FaceAlpha',1., 'LineWidth',4.);

r = 0
saveas(gcf,['./Figures/3d_field_idealized' num2str(m) '__id_' num2str(r) '_.png']);

% -------------------------------------------------------------
% -------------------------------------------------------------

function d = can_distance(X, Y)
    global Rhom
    Z = mod(X - Y,1); 
    Z = Rhom(:,[1; 3])*Z;
    n1 = vecnorm(Z - Rhom(:,1), 2, 1);
    n2 = vecnorm(Z - Rhom(:,2), 2, 1);
    n3 = vecnorm(Z - Rhom(:,3), 2, 1);
    n4 = vecnorm(Z - Rhom(:,4), 2, 1);
    D = vertcat(n1,n2,n3,n4);
    d = min(D,[],1);
end


function g = normalize(f)
    g_ = f - min(f,[],'all');
    g_ = g_./max(g_, [],'all');
    g  = g_;
end

function M=RandOrthMat(n, tol)
% M = RANDORTHMAT(n)
% generates a random n x n orthogonal real matrix.
%
% M = RANDORTHMAT(n,tol)
% explicitly specifies a thresh value that measures linear dependence
% of a newly formed column with the existing columns. Defaults to 1e-6.
%
% In this version the generated matrix distribution *is* uniform over the manifold
% O(n) w.r.t. the induced R^(n^2) Lebesgue measure, at a slight computational 
% overhead (randn + normalization, as opposed to rand ). 
% 
% (c) Ofek Shilon , 2006.
    if nargin==1
	  tol=1e-6;
    end
    
    M = zeros(n); % prealloc
    
    % gram-schmidt on random column vectors
    
    vi = randn(n,1);  

    % the n-dimensional normal distribution has spherical symmetry, which implies
    % that after normalization the drawn vectors would be uniformly distributed on the
    % n-dimensional unit sphere.
    M(:,1) = vi ./ norm(vi);
    
    for i=2:n
	  nrm = 0;
	  while nrm<tol
		vi = randn(n,1);
		vi = vi -  M(:,1:i-1)  * ( M(:,1:i-1).' * vi )  ;
		nrm = norm(vi);
	  end
	  M(:,i) = vi ./ nrm;
    end %i
        
end  % RandOrthMat


function g=get_bubbles(n, X)

    global Rhom
    Rhom = [cos(0) sin(0); cos(pi/3.)+cos(0) sin(pi/3.)+sin(0); cos(pi/3.) sin(pi/3.); 0 0]';

    
    A = zeros(n, 2, 3);
    B = zeros(n, 3, 3);
    F = zeros(n, size(X,2));
    for i = 1:n
    % ---------------------
    B(i,[1;2],[1;2]) = Rhom(:,[1; 3]);
    end


    for i = 1:n
    % ---------------------

    B(i,:,3) = [0;0;1];
    B = B.*(.9 + (rand-0.5)*2*0.1);
    Ai = inv(reshape(B(i,1:3,1:3),3,3));
    A(i,:,:) = Ai(1:2,:)*RandOrthMat(3);
    end

    % A(1,:,:) = A1;
    % A(2,:,:) = A6;
    % A(3,:,:) = A3;
    % A(4,:,:) = A4;
    % A(3,:,:) = A5;
    % A(4,:,:) = A3;

    for i = 1:n
    Ai  = reshape(A(i,:,:), 2,3);
    fi  = can_distance(Ai*X, zeros(2,size(X,2)) );
    fi   = exp(-fi.^2/1.);
    F(i,:) = fi;

    % ---------------------
    end
    

    f = F(1,:);
    for i = 1:n
        f  = f + F(i,:);
    end
    
    g = f;
    g = g/sum(g);
        
end




