%%%%tests prélimimaires de validation des algos de krigeage
%%avec élimination de voxels connus puis leur interpolation à posteriori


%% tests à l'echelle de la branche avec voxels de 18cm de coté

load('data_branch_2krig.mat')

%%Experience 281, Br1, 2.5m, Feuillée
nbVox=length(indicesSerreLVOX_varBulk.N{281})%64
base=nbVox^(1/3)
reso=0.7/base

SEQ= reshape([reso+reso*0.5:reso:reso*base+reso*0.5],length(reso+reso*0.5:reso:reso*base+reso*0.5),1);

%creation du système de coordonnées selon ordre indices
RDI_coordX=repmat(SEQ,floor(nbVox/base),1);
RDI_coordY=repmat(SEQ,floor(nbVox/base^2),floor(nbVox/base^2))';RDI_coordY=RDI_coordY(:)';RDI_coordY=reshape(RDI_coordY,length(RDI_coordY),1);
RDI_coordZ=repmat(SEQ,1,floor(nbVox/base))';RDI_coordZ=RDI_coordZ(:)';RDI_coordZ=reshape(RDI_coordZ,length(RDI_coordZ),1);
RDI_coord=[RDI_coordX RDI_coordY RDI_coordZ];clear RDI_coordX RDI_coordY RDI_coordZ


RDI=reshape(indicesSerreLVOX_varBulk.RDI{281},length(indicesSerreLVOX_varBulk.RDI{281}),1);%import des RDI

RDI_all=[RDI_coord RDI];clear RDI

%REmove voxels outside sphere r=0.35 for experimental variogram computation
ind=find(sqrt((repmat(0.35,nbVox,1)-RDI_all(:,1)).^2+(repmat(0.35,nbVox,1)-RDI_all(:,2)).^2+(repmat(0.35,nbVox,1)-RDI_all(:,3)).^2)<=0.35);
RDI=RDI_all(ind,:);
%REmove Nan RDI
ind=find(~isnan(RDI(:,4)));
RDI=RDI(ind,:); clear ind

%Remove choosen voxel(s) : n=1 /!\  voxel n=1 at sphere limit (might be
%overestimated)
%n=20, moyennement plein
%n=13, très plein
n=20
Exp_I_k=RDI(n,4);
coord_excl=RDI(n,1:3);
RDI=RDI([1:n-1,n+1:20],:);


%Experimental Variogram computing
RDIvgm=variogram(RDI(:,1:3),RDI(:,4),'nrbins',5,'maxdist',0.35,'plotit',true)

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


%Computing distance matrix for measured RDI
distsObs=distmat(RDI(:,1:3));
%Computing covariance matrix for measured RDI
COVij=(RDIfit.sill-RDIfit.nugget).*(1-exp(-distsObs./(RDIfit.range)));
%Add terms to Cij matrix for sum(alpha)=1
COVij=[COVij repmat(1,length(COVij),1)];
COVij=[COVij;[repmat(1,1,(length(COVij)-1)) 0]];

%%% POUR UN SEUL POINT A KRIGER (ORDINAIRE)

%Computing distance matrix between measured and one unknown
dists0i_1=distmat(RDI(:,1:3),coord_excl(1,:))
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

%Evaluation du krigeage
I_k/Exp_I_k% n=1 très surestimé (comme prévu car bordure)
%n=20 : 21% de surestimation
%n=13: 58% de sousestimation
%A noter que ces résultats moyens (mais du même ordre de grandeur) peuvent
%être expliqués d'une part à cause du faible nombre de valeurs pour le fit
%et à cause de la forte heterogeneite spatiale (par exemple le voxel n=13
%très plein et entourré de voxel nettement plus vides


%% tests à l'echelle de l'arbre avec voxels de ??