x = rand(1000,1)*4-2;  
y = rand(1000,1)*4-2;
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
d = variogram([x y],z,'plotit',true,'nrbins',50);
title('Isotropic variogram')
subplot(2,2,4)
d2 = variogram([x y],z,'plotit',true,'nrbins',10,'anisotropy',true);
title('Anisotropic variogram')

vgram = d2;

a0 = 0.1;
c0 = 5;
fit = struct();
[fit.range fit.sill fit.nugget] = variogramfit(vgram.distance, vgram.val, a0, c0, [], 'plotit', false);

vx = 0:0.05:2;
vy = (fit.sill).*(1-exp(-vx./(fit.range)));
figure
hold on;
plot(vx, vy);
plot(vgram.distance, vgram.val);
hold off

distsObs = distmat(x);

COVij=(fit.sill).*(1-exp(-distsObs./(fit.range)));

% la matrice des distances n'est pas triée
COVij=[COVij repmat(1,length(COVij),1)];
COVij=[COVij;[repmat(1,1,(length(COVij)-1)) 0]];

ux = 1.5;
%ux = (pi/(2*15));
uy = ux;

z_exp = 3*sin(ux*15);
%z_exp = 3*ux + uy;

dist_uxy = distmat([x y], [ux uy]);

cov=(fit.sill).*(1-exp(-dist_uxy./(fit.range*3)));

cov=[cov;1];

w = inv(COVij) * cov;

w=w(1:length(w)-1,:);

z_act=w'*z;

error = abs(z_act - z_exp)
