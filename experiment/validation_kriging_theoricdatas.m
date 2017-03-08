%%% test de l'algo de kriging sur un jeu de donnée théorique controlé

x = rand(1000,1)*4-2;  
y = rand(1000,1)*4-2;
coord=[x y];
z = 3*sin(x*15) + randn(size(x));
subplot(2,2,1)
scatter(x,y,4,z,'filled'); box on;
ylabel('y'); xlabel('x')
title('data (coloring according to z-value)')
subplot(2,2,2)
hist(z,20)
ylabel('frequency'); xlabel('z')
title('histogram of z-values')
subplot(2,2,3)
d = variogram(coord,z,'plotit',true,'nrbins',50);
title('Isotropic variogram')
subplot(2,2,4)
d2 = variogram([x y],z,'plotit',true,'nrbins',10,'anisotropy',true);
title('Anisotropic variogram')

vgram = d2;

a0 = 0.1;
c0 = 5;
fit = struct();
%Experimental Variogram computing
figure();
d = variogram(coord,z,'plotit',true,'nrbins',50);

%Fitting theoretical variogram
hold on
h=vgram.distance;
gammaexp=vgram.val;
numobs=vgram.num;
fit=struct();
[fit.range fit.sill fit.nugget fit.S]=variogramfit(h,gammaexp,a0,c0,numobs,'model','exponential','nugget',0.00,'weightfun','none','plotit',true)
clear gammaexp h numobs a0 c0
%[fit.range fit.sill fit.nugget] = variogramfit(vgram.distance, vgram.val, a0, c0, [], 'plotit', true);


vx = 0:0.05:2;
vy = (fit.sill).*(1-exp(-vx./(fit.range)));
figure
hold on;
plot(vx, vy);
plot(vgram.distance, vgram.val);
hold off

%Computing distance matrix for measured x
distsObs = distmat(coord);% la matrice des distances n'est pas triée

COVij=(fit.sill).*(1-exp(-distsObs./(fit.range)));
%Add terms to Cij matrix for sum(alpha)=1
COVij=[COVij repmat(1,length(COVij),1)];
COVij=[COVij;[repmat(1,1,(length(COVij)-1)) 0]];

%%% POUR UN SEUL POINT A KRIGER (ORDINAIRE)
ux = (pi/(2*15));
uy = ux;
coord_2krig=[ux uy];

z_exp = 3*sin(ux*15)

%Computing distance matrix between measured and one unknown
dist_uxy = distmat(coord, coord_2krig);

%Computing covariance matrix between measured and one unknown
cov=(fit.sill).*(1-exp(-dist_uxy./(fit.range)));

cov=[cov;1];

w = inv(COVij) * cov;

w=w(1:length(w)-1,:);
sum(w)%ok

z_act=w'*z;

error = abs(z_act - z_exp)