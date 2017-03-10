%%% test experimental variogram computing


load('data_branch_2krig.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ORDINARY KRIGING

%create input datas
nbVox=4096;
base=nbVox^(1/3)
reso=0.7/base;%sphere de largeur 0.7


SEQ= reshape([reso+reso*0.5:reso:reso*base+reso*0.5],length(reso+reso*0.5:reso:reso*base+reso*0.5),1);

%creation du système de coordonnées selon ordre indices
RDI_coordX=repmat(SEQ,floor(nbVox/base),1);
RDI_coordY=repmat(SEQ,floor(nbVox/base^2),floor(nbVox/base^2))';RDI_coordY=RDI_coordY(:)';RDI_coordY=reshape(RDI_coordY,length(RDI_coordY),1);
RDI_coordZ=repmat(SEQ,1,floor(nbVox/base))';RDI_coordZ=RDI_coordZ(:)';RDI_coordZ=reshape(RDI_coordZ,length(RDI_coordZ),1);
RDI_coord=[RDI_coordX RDI_coordY RDI_coordZ];clear RDI_coordX RDI_coordY RDI_coordZ

RDI=reshape(indicesSerreLVOX_varBulk.RDI{1},length(indicesSerreLVOX_varBulk.RDI{1}),1);%import des RDI

RDI_all=[RDI_coord RDI];clear RDI


%REmove Nan RDI (OCCLUDED, TO KRIGE)
ind=find(~isnan(RDI_all(:,4)));
RDI=RDI_all(ind,:); %582 voxels occluded in all the cubic scene
%outside sphere:1293<0.01 + 158 entre 0.01 et 0.6, dont 96<0.1 + 539=NaN (-->occlus)

%REmove voxels outside sphere r=0.35 (TO UNBIAS VARIOG COMPUTATION)
ind=find(sqrt((repmat(0.35,length(RDI),1)-RDI(:,1)).^2+(repmat(0.35,length(RDI),1)-RDI(:,2)).^2+(repmat(0.35,length(RDI),1)-RDI(:,3)).^2)<=0.35);
RDI=RDI(ind,:);%42 voxels occlus dans la sphere

%Coord of occluded voxel
ind2=find(isnan(RDI_all(:,4)));
RDI_2krig=RDI_all(ind2,:);%ceux que l'on va kriger

%Coord of voxel with measured RDI outside of vegetation sphere
ind1=find(sqrt((repmat(0.35,nbVox,1)-RDI_all(:,1)).^2+(repmat(0.35,nbVox,1)-RDI_all(:,2)).^2+(repmat(0.35,nbVox,1)-RDI_all(:,3)).^2)>0.35);
ind2=find(~isnan(RDI_all(ind1,4)));
RDI_outsideObs=RDI_all(ind1(ind2),:);


clear ind1 ind2



%Experimental Variogram computing
RDIvgm=variogram(RDI(:,1:3),RDI(:,4),'nrbins',15,'maxdist',0.35,'plotit',true)

%Fitting theoretical variogram
hold on
a0=0.2;%range
c0=0.018;%sill
h=RDIvgm.distance;
gammaexp=RDIvgm.val;
numobs=RDIvgm.num;
RDIfit=struct();
[RDIfit.range RDIfit.sill RDIfit.nugget RDIfit.S]=variogramfit(h,gammaexp,a0,c0,numobs,'model','exponential','nugget',0.02,'weightfun','none','plotit',true)
clear gammaexp h numobs a0 c0

%%%Computing covariance matrix from fitted empirical variogram

%Ajout des RDI nuls observés en dehors de la sphere
RDI=[RDI;RDI_outsideObs];

%Computing distance matrix for measured RDI
distsObs=distmat(RDI(:,1:3));
%Computing covariance matrix for measured RDI
COVij=(RDIfit.sill-RDIfit.nugget).*(1-exp(-distsObs./(RDIfit.range)));
%Add terms to Cij matrix for sum(alpha)=1
COVij=[COVij repmat(1,length(COVij),1)];
COVij=[COVij;[repmat(1,1,(length(COVij)-1)) 0]];


%%% POUR UN SEUL POINT A KRIGER (ORDINAIRE)

%Computing distance matrix between measured and one unknown
dists0i_1=distmat(RDI(:,1:3),RDI_2krig(1,1:3))
%Computing covariance matrix between measured and one unknown
COV0i_1=(RDIfit.sill-RDIfit.nugget).*(1-exp(-dists0i_1./(RDIfit.range)));
%Add terms 1 to C0i matrix for sum(alpha)=1
COV0i_1=[COV0i_1;1];

%COmputing weights matrix (alpha)
Weights=inv(COVij)*COV0i_1;
%Remove last weight (artificial term)
Weights=Weights(1:length(Weights)-1,:);
sum(Weights)%check =1

%Compute Weights transpose
Weights=(Weights)';

%Computing ~I at point RDI_unk_coord(1,:)
I_k=Weights*RDI(:,4)

%%% POUR PLUSIEURS POINT A KRIGER (ORDINAIRE)

%Computing distance matrix between measured and unknown
dists0i_2krig=[];
for i=1:length(RDI_2krig)
dists0i_2krig(:,i)=distmat(RDI(:,1:3),RDI_2krig(i,1:3));
end

%Computing covariance matrix between measured and one unknown
COV0i=[];
for i=1:length(RDI_2krig)
COV0i(:,i)=(RDIfit.sill-RDIfit.nugget).*(1-exp(-dists0i_2krig(:,i)./(RDIfit.range)));
end

%Add terms 1 to C0i matrix for sum(alpha)=1
COV0i=[COV0i;repmat(1,1,length(RDI_2krig))];

%COmputing weights matrix (alpha)
Weights=inv(COVij)*COV0i;
%Remove last weight (artificial term)
Weights=Weights(1:length(Weights)-1,:);
unique(sum(Weights))%check : ok

%Compute Weights transpose
Weights=(Weights)';

%Computing ~I at points RDI_unk_coord(:,:)
I_k=Weights*RDI(:,4);
%REplace negative values by 0
ind=find(I_k<0)
I_k(ind)=0;

RDI_kriged=[RDI_2krig(:,1:3), I_k];

%% Visualisation mesurés/krigés ordinaires

% [X1 Y1]=meshgrid(RDI_all(ind,1),RDI_all(ind,2));
% Z1=griddata(RDI_all(ind,1),RDI_all(ind,2),RDI_all(ind,4), X1,Y1);
% figure();surface(X1, Y1, Z1)
figure();
ind= find(RDI_kriged(:,3)==0.328125);
scatter(RDI_kriged(ind,1),RDI_kriged(ind,2),(RDI_kriged(ind,4).*100+0.000000001),'b')
xlim([0,0.8]);ylim([0,0.8]);

hold on
ind= find(RDI(:,3)==0.328125);
scatter(RDI(ind,1),RDI(ind,2),(RDI(ind,4).*100+0.000000001),'r','+')
title('blue= ordinary kriged red=measured');xlim([0,0.8]);ylim([0,0.8]);
hold off 
clear ind

figure();
scatter3(RDI_kriged(:,1),RDI_kriged(:,2),RDI_kriged(:,3),(RDI_kriged(:,4).*100+0.000000001),'b')
xlim([0,0.8]);ylim([0,0.8]);zlim([0,0.8])
hold on
scatter3(RDI(:,1),RDI(:,2),RDI(:,3),(RDI(:,4).*100+0.000000001),'r')
title('blue= ordinary kriged red=measured');xlim([0,0.8]);ylim([0,0.8]);zlim([0,0.8])
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BINOMIAL KRIGING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% BINOMIAL WITH ORDINARY EXPERIMENTAL VARIOGRAM
%create input datas
nbVox=4096;
base=nbVox^(1/3)
reso=0.7/base;%sphere de largeur 0.7


SEQ= reshape([reso+reso*0.5:reso:reso*base+reso*0.5],length(reso+reso*0.5:reso:reso*base+reso*0.5),1);

%creation du système de coordonnées selon ordre indices
RDI_coordX=repmat(SEQ,floor(nbVox/base),1);
RDI_coordY=repmat(SEQ,floor(nbVox/base^2),floor(nbVox/base^2))';RDI_coordY=RDI_coordY(:)';RDI_coordY=reshape(RDI_coordY,length(RDI_coordY),1);
RDI_coordZ=repmat(SEQ,1,floor(nbVox/base))';RDI_coordZ=RDI_coordZ(:)';RDI_coordZ=reshape(RDI_coordZ,length(RDI_coordZ),1);
RDI_coord=[RDI_coordX RDI_coordY RDI_coordZ];clear RDI_coordX RDI_coordY RDI_coordZ

RDI=reshape(indicesSerreLVOX_varBulk.RDI{1},length(indicesSerreLVOX_varBulk.RDI{1}),1);%import des RDI
RDI_N=reshape(indicesSerreLVOX_varBulk.N{1},length(indicesSerreLVOX_varBulk.N{1}),1);%import des Nt-Nb

RDI_all=[RDI_coord RDI RDI_N];clear RDI RDI_N


%REmove Nan RDI (OCCLUDED, TO KRIGE)
ind=find(~isnan(RDI_all(:,4)));
RDI=RDI_all(ind,:); %582 voxels occluded in all the cubic scene
%outside sphere:1293<0.01 + 158 entre 0.01 et 0.6, dont 96<0.1 + 539=NaN (-->occlus)

%REmove voxels outside sphere r=0.35 (TO UNBIAS VARIOG COMPUTATION)
ind=find(sqrt((repmat(0.35,length(RDI),1)-RDI(:,1)).^2+(repmat(0.35,length(RDI),1)-RDI(:,2)).^2+(repmat(0.35,length(RDI),1)-RDI(:,3)).^2)<=0.35);
RDI=RDI(ind,:);%42 voxels occlus dans la sphere

%Coord of occluded voxel
ind2=find(isnan(RDI_all(:,4)));
RDI_2krig=RDI_all(ind2,:);%ceux que l'on va kriger

%Coord of voxel with measured RDI outside of vegetation sphere
ind1=find(sqrt((repmat(0.35,nbVox,1)-RDI_all(:,1)).^2+(repmat(0.35,nbVox,1)-RDI_all(:,2)).^2+(repmat(0.35,nbVox,1)-RDI_all(:,3)).^2)>0.35);
ind2=find(~isnan(RDI_all(ind1,4)));
RDI_outsideObs=RDI_all(ind1(ind2),:);


clear ind1 ind2



%Experimental Variogram computing NON BINO
RDIvgm=variogram(RDI(:,1:3),RDI(:,4),'nrbins',15,'maxdist',0.35,'plotit',true)

%Fitting theoretical variogram
hold on
a0=0.2;%range
c0=0.018;%sill
h=RDIvgm.distance;
gammaexp=RDIvgm.val;
numobs=RDIvgm.num;
RDIfit=struct();
[RDIfit.range RDIfit.sill RDIfit.nugget RDIfit.S]=variogramfit(h,gammaexp,a0,c0,numobs,'model','exponential','nugget',0.02,'weightfun','none','plotit',true)
clear gammaexp h numobs a0 c0

%%%Computing covariance matrix from fitted empirical variogram

%Ajout des RDI nuls observés en dehors de la sphere
RDI=[RDI;RDI_outsideObs];

%Computing distance matrix for measured RDI
distsObs=distmat(RDI(:,1:3));
%Computing covariance matrix for measured RDI
COVij=(RDIfit.sill-RDIfit.nugget).*(1-exp(-distsObs./(RDIfit.range)));
%Compute mean of I_chap
ind=find(RDI(:,5)>30);pi_hat=mean(RDI(ind,4));
%Compute variance of I_chap
ind=find(RDI(:,5)>30);sig2_hat=var(RDI(ind,4));
%Add BINOMIAL term to diagonal of COVij TO MODIFY alphai
COVij=COVij+diag((1./RDI(:,5)).*(pi_hat*(1-pi_hat)-sig2_hat));

%Add terms to Cij matrix for sum(alpha)=1
COVij=[COVij repmat(1,length(COVij),1)];
COVij=[COVij;[repmat(1,1,(length(COVij)-1)) 0]];


%%% POUR UN SEUL POINT A KRIGER (ORDINAIRE)

%Computing distance matrix between measured and one unknown
dists0i_1=distmat(RDI(:,1:3),RDI_2krig(1,1:3))
%Computing covariance matrix between measured and one unknown
COV0i_1=(RDIfit.sill-RDIfit.nugget).*(1-exp(-dists0i_1./(RDIfit.range)));
%Add terms 1 to C0i matrix for sum(alpha)=1
COV0i_1=[COV0i_1;1];

%COmputing weights matrix (alpha)
Weights=inv(COVij)*COV0i_1;
%Remove last weight (artificial term)
Weights=Weights(1:length(Weights)-1,:);
sum(Weights)%check =1

%Compute Weights transpose
Weights=(Weights)';

%Computing ~I at point RDI_unk_coord(1,:)
I_k_bino_novgm=Weights*RDI(:,4)

%%% POUR PLUSIEURS POINT A KRIGER (BINOMIAL WITHOUT VARIOGRAM)

%Computing distance matrix between measured and unknown
dists0i_2krig=[];
for i=1:length(RDI_2krig)
dists0i_2krig(:,i)=distmat(RDI(:,1:3),RDI_2krig(i,1:3));
end

%Computing covariance matrix between measured and one unknown
COV0i=[];
for i=1:length(RDI_2krig)
COV0i(:,i)=(RDIfit.sill-RDIfit.nugget).*(1-exp(-dists0i_2krig(:,i)./(RDIfit.range)));
end

%Add terms 1 to C0i matrix for sum(alpha)=1
COV0i=[COV0i;repmat(1,1,length(RDI_2krig))];

%COmputing weights matrix (alpha)
Weights=inv(COVij)*COV0i;
%Remove last weight (artificial term)
Weights=Weights(1:length(Weights)-1,:);
unique(sum(Weights))%check : ok

%Compute Weights transpose
Weights=(Weights)';

%Computing ~I at points RDI_unk_coord(:,:)
I_k_bino_novgm=Weights*RDI(:,4);
%REplace negative values by 0
ind=find(I_k_bino_novgm<0)
I_k_bino_novgm(ind)=0;

RDI_kriged_bino_novgm=[RDI_2krig(:,1:3), I_k];



%% BINOMIAL WITH BINOMIAL EXPERIMENTAL VARIOGRAM (FULL BINOMIAL METHOD)

%create input datas
nbVox=4096;
base=nbVox^(1/3)
reso=0.7/base;%sphere de largeur 0.7


SEQ= reshape([reso+reso*0.5:reso:reso*base+reso*0.5],length(reso+reso*0.5:reso:reso*base+reso*0.5),1);

%creation du système de coordonnées selon ordre indices
RDI_coordX=repmat(SEQ,floor(nbVox/base),1);
RDI_coordY=repmat(SEQ,floor(nbVox/base^2),floor(nbVox/base^2))';RDI_coordY=RDI_coordY(:)';RDI_coordY=reshape(RDI_coordY,length(RDI_coordY),1);
RDI_coordZ=repmat(SEQ,1,floor(nbVox/base))';RDI_coordZ=RDI_coordZ(:)';RDI_coordZ=reshape(RDI_coordZ,length(RDI_coordZ),1);
RDI_coord=[RDI_coordX RDI_coordY RDI_coordZ];clear RDI_coordX RDI_coordY RDI_coordZ

RDI=reshape(indicesSerreLVOX_varBulk.RDI{1},length(indicesSerreLVOX_varBulk.RDI{1}),1);%import des RDI
RDI_N=reshape(indicesSerreLVOX_varBulk.N{1},length(indicesSerreLVOX_varBulk.N{1}),1);%import des Nt-Nb

RDI_all=[RDI_coord RDI RDI_N];clear RDI RDI_N


%REmove Nan RDI (OCCLUDED, TO KRIGE)
ind=find(~isnan(RDI_all(:,4)));
RDI=RDI_all(ind,:); %582 voxels occluded in all the cubic scene
%outside sphere:1293<0.01 + 158 entre 0.01 et 0.6, dont 96<0.1 + 539=NaN (-->occlus)

%REmove voxels outside sphere r=0.35 (TO UNBIAS VARIOG COMPUTATION)
ind=find(sqrt((repmat(0.35,length(RDI),1)-RDI(:,1)).^2+(repmat(0.35,length(RDI),1)-RDI(:,2)).^2+(repmat(0.35,length(RDI),1)-RDI(:,3)).^2)<=0.35);
RDI=RDI(ind,:);%42 voxels occlus dans la sphere

%Coord of occluded voxel
ind2=find(isnan(RDI_all(:,4)));
RDI_2krig=RDI_all(ind2,:);%ceux que l'on va kriger

%Coord of voxel with measured RDI outside of vegetation sphere
ind1=find(sqrt((repmat(0.35,nbVox,1)-RDI_all(:,1)).^2+(repmat(0.35,nbVox,1)-RDI_all(:,2)).^2+(repmat(0.35,nbVox,1)-RDI_all(:,3)).^2)>0.35);
ind2=find(~isnan(RDI_all(ind1,4)));
RDI_outsideObs=RDI_all(ind1(ind2),:);

clear ind1 ind2

%Experimental BINOMIAL Variogram computing
RDIvgm=variogramVOX_beta(RDI(:,1:3),RDI(:,4),RDI(:,5),'nrbins',15,'maxdist',0.35,'binomial',1,'plotit',true)

%Fitting theoretical variogram
hold on
a0=0.2;%range
c0=0.018;%sill
h=RDIvgm.distance;
gammaexp=RDIvgm.val;
numobs=RDIvgm.num;
RDIfit=struct();
[RDIfit.range RDIfit.sill RDIfit.nugget RDIfit.S]=variogramfit(h,gammaexp,a0,c0,numobs,'model','exponential','nugget',0.02,'weightfun','none','plotit',true)
clear gammaexp h numobs a0 c0

%%%Computing covariance matrix from fitted empirical variogram

%Ajout des RDI nuls observés en dehors de la sphere
RDI=[RDI;RDI_outsideObs];

%Computing distance matrix for measured RDI
distsObs=distmat(RDI(:,1:3));
%Computing covariance matrix for measured RDI
COVij=(RDIfit.sill-RDIfit.nugget).*(1-exp(-distsObs./(RDIfit.range)));
%Compute mean of I_chap
ind=find(RDI(:,5)>30);pi_hat=mean(RDI(ind,4));
%Compute variance of I_chap
ind=find(RDI(:,5)>30);sig2_hat=var(RDI(ind,4));
%Add BINOMIAL term to diagonal of COVij TO MODIFY alphai
COVij=COVij+diag((1./RDI(:,5)).*(pi_hat*(1-pi_hat)-sig2_hat));

%Add terms to Cij matrix for sum(alpha)=1
COVij=[COVij repmat(1,length(COVij),1)];
COVij=[COVij;[repmat(1,1,(length(COVij)-1)) 0]];


%%% POUR UN SEUL POINT A KRIGER (FULL BINOMIAL METHOD)

%Computing distance matrix between measured and one unknown
dists0i_1=distmat(RDI(:,1:3),RDI_2krig(1,1:3))
%Computing covariance matrix between measured and one unknown
COV0i_1=(RDIfit.sill-RDIfit.nugget).*(1-exp(-dists0i_1./(RDIfit.range)));
%Add terms 1 to C0i matrix for sum(alpha)=1
COV0i_1=[COV0i_1;1];

%COmputing weights matrix (alpha)
Weights=inv(COVij)*COV0i_1;
%Remove last weight (artificial term)
Weights=Weights(1:length(Weights)-1,:);
sum(Weights)%check =1

%Compute Weights transpose
Weights=(Weights)';

%Computing ~I at point RDI_unk_coord(1,:)
I_k_bino_novgm=Weights*RDI(:,4)

%%% POUR PLUSIEURS POINT A KRIGER (ORDINAIRE)

%Computing distance matrix between measured and unknown
dists0i_2krig=[];
for i=1:length(RDI_2krig)
dists0i_2krig(:,i)=distmat(RDI(:,1:3),RDI_2krig(i,1:3));
end

%Computing covariance matrix between measured and one unknown
COV0i=[];
for i=1:length(RDI_2krig)
COV0i(:,i)=(RDIfit.sill-RDIfit.nugget).*(1-exp(-dists0i_2krig(:,i)./(RDIfit.range)));
end

%Add terms 1 to C0i matrix for sum(alpha)=1
COV0i=[COV0i;repmat(1,1,length(RDI_2krig))];

%COmputing weights matrix (alpha)
Weights=inv(COVij)*COV0i;
%Remove last weight (artificial term)
Weights=Weights(1:length(Weights)-1,:);
unique(sum(Weights))%check : ok

%Compute Weights transpose
Weights=(Weights)';

%Computing ~I at points RDI_unk_coord(:,:)
I_k_bino=Weights*RDI(:,4);
%REplace negative values by 0
ind=find(I_k_bino<0)
I_k_bino(ind)=0;

RDI_kriged_bino=[RDI_2krig(:,1:3), I_k];

%% Visualisation mesurés/krigés binomial

% [X1 Y1]=meshgrid(RDI_all(ind,1),RDI_all(ind,2));
% Z1=griddata(RDI_all(ind,1),RDI_all(ind,2),RDI_all(ind,4), X1,Y1);
% figure();surface(X1, Y1, Z1)
figure();
ind= find(RDI_kriged_bino(:,3)==0.328125);
scatter(RDI_kriged_bino(ind,1),RDI_kriged_bino(ind,2),(RDI_kriged_bino(ind,4).*100+0.000000001),'b')
xlim([0,0.8]);ylim([0,0.8]);

hold on
ind= find(RDI(:,3)==0.328125);
scatter(RDI(ind,1),RDI(ind,2),(RDI(ind,4).*100+0.000000001),'r','+')
title('blue= ordinary kriged red=measured');xlim([0,0.8]);ylim([0,0.8]);
hold off 
clear ind

figure();
scatter3(RDI_kriged_bino(:,1),RDI_kriged_bino(:,2),RDI_kriged_bino(:,3),(RDI_kriged_bino(:,4).*100+0.000000001),'b')
xlim([0,0.8]);ylim([0,0.8]);zlim([0,0.8])
hold on
scatter3(RDI(:,1),RDI(:,2),RDI(:,3),(RDI(:,4).*100+0.000000001),'r')
title('blue= ordinary kriged red=measured');xlim([0,0.8]);ylim([0,0.8]);zlim([0,0.8])
hold off

%Comparaison ordinaire/krigé
figure();plot(I_k,I_k_bino,'+')
sum(I_k)
sum(I_k_bino_novgm)
sum(I_k_bino)