%%% Adaptation de la fonction variogram en vue d'un krigeage binomial
%9/03/2017
% Max, from Wolfgang Schwanghart

%% loading test dataset
load('data_branch_2krig.mat')

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

%REmove Nan RDI (OCCLUDED)
ind=find(~isnan(RDI_all(:,4)));
RDI=RDI_all(ind,:); %582 voxels occluded in all the cubic scene
%outside sphere:1293<0.01 + 158 entre 0.01 et 0.6, dont 96<0.1 + 539=NaN (-->occlus)

%REmove voxels outside sphere r=0.35 (TO UNBIAS VARIOG COMPUTATION)
ind=find(sqrt((repmat(0.35,length(RDI),1)-RDI(:,1)).^2+(repmat(0.35,length(RDI),1)-RDI(:,2)).^2+(repmat(0.35,length(RDI),1)-RDI(:,3)).^2)<=0.35);
RDI=RDI(ind,:);%42 voxels occlus dans la sphere



%% Experimental variogram computing
clear S iid

%% Arguments settings
x=[RDI(:,1),RDI(:,2),RDI(:,3)];%matrice de coordonnées des points observés
y=RDI(:,4);%variable RDI
N=RDI(:,5);%variable N=Nt-Nb

nrbins=15;
maxdist=0.35;
binomial=1;
plotit=true;

varargin={'nrbins',nrbins,'maxdist',maxdist,'binomial',binomial,'plotit',plotit};%list of non-default arguments values


%% size checking
if size(y,1) ~= size(x,1);
    error('x and y must have the same number of rows')
end

%% check for nans
II = any(isnan(x),2) | isnan(y);%tests matrixs rows
x(II,:) = [];
y(II)   = [];
if unique(II) ~= 0;
warning('NaN values have been removed')
end

%% extent of dataset
minx = min(x,[],1);
maxx = max(x,[],1);
maxd = sqrt(sum((maxx-minx).^2));
nrdims = size(x,2);

%% check input using PARSEARGS (parameters values)
%default values:
params.nrbins      = 20;
params.maxdist     = maxd/2;
params.type        = {'default','gamma','cloud1','cloud2'};
params.plotit      = false;
params.anisotropy  = false;
params.thetastep   = 30;
params.subsample   = inf;
params.binomial    = 0;

params = parseargs(params,varargin{:});%replace default values if specified

if params.maxdist > maxd;%check specified maximum dist
    warning('Matlab:Variogram',...
            ['Maximum distance exceeds maximum distance \n' ... 
             'in the dataset. maxdist was decreased to ' num2str(maxd) ]);
    params.maxdist  = maxd;
end

if params.anisotropy && nrdims ~= 2 %%% TO ADAPT
    params.anisotropy = false;
    warning('Matlab:Variogram',...
            'Anistropy is only supported for 2D data');
end

%% if specified, take only a subset of the data at random;
if ~isinf(params.subsample) && numel(y)>params.subsample;
    IX = randperm(numel(y),params.subsample);
    x  = x(IX,:);
    y  = y(IX,:);
end

%% calculate bin tolerance
tol = params.maxdist/params.nrbins; %pas de taille des classes

%% calculate constrained distance matrix
iid = distmatBound(x,params.maxdist);

%% calculate squared difference between values of coordinate pairs BINOMIAL MODIFIED
if params.binomial
    Npair_product=RDI(iid(:,1),5).*RDI(iid(:,2),5);
    Npair_sum=RDI(iid(:,1),5)+RDI(iid(:,2),5);
    ind=find(RDI(:,5)>=30);%TO INCREASE TO BE MORE CONSERVATIVE
    expct=mean(RDI(ind,4));%OBSERVED expectation of RDI
    var2=var(RDI(ind,4));%OBSERVED variance of RDI
    Nterm=Npair_product./Npair_sum;
    
        lam = (Nterm.*((y(iid(:,1))-y(iid(:,2))).^2)-expct*(1-expct)+var2);%/!\ valeurs légèrement négative
else
        lam      = (y(iid(:,1))-y(iid(:,2))).^2;
end
    iid=[iid lam];

%% calculate variogram
switch params.type
    case {'default','gamma'}
            % variogram anonymous function
        fvar     = @(x) 1./(2*numel(x)) * sum(x);%s'applique à lam
        
        % distance bins
        edges      = linspace(0,params.maxdist,params.nrbins+1);%distance classes
        edges(end) = inf;

        [nedge,ixedge] = histc(iid(:,3),edges);%nedge=nb of pairs by class
        %ixedge=attributed class
        iid=[iid ixedge];
        
        if params.anisotropy
            S.val      = accumarray([ixedge ixtheta],lam,...
                                 [numel(edges) numel(thetaedges)],fvar,nan);
            S.val(:,end)=S.val(:,1); 
            S.theta    = thetacents;
            S.num      = accumarray([ixedge ixtheta],ones(size(lam)),...
                                 [numel(edges) numel(thetaedges)],@sum,nan);
            S.num(:,end)=S.num(:,1);                 
        else
            if params.binomial
               for i=2:params.nrbins
                   ind=find(iid(:,5)==i);
                  S.W(i)=sum(Npair_product(ind)./Npair_sum(ind));
               end 
               S.W(1)=NaN;
               S.W(end+1)=NaN;
               S.W=(S.W)';
               S.val=accumarray(ixedge,lam,[numel(edges) 1],@sum,nan);
               S.val=(1./(2.*S.W)).*S.val;
            else
            S.val      = accumarray(ixedge,lam,[numel(edges) 1],fvar,nan);%compute gamma by distance class  
            %apply fvar to the values in lam that have identical subscripts in ixedge
            %output of size numel(edges)
            end
            S.num      = accumarray(ixedge,ones(size(lam)),[numel(edges) 1],@sum,nan);%count pairs nb
        end
        S.distance = (edges(1:end-1)+tol/2)';
        S.val(end,:) = [];
        S.num(end,:) = [];

    case 'cloud1'
        edges      = linspace(0,params.maxdist,params.nrbins+1);
        edges(end) = inf;
        
        [nedge,ixedge] = histc(iid(:,3),edges);
        
        S.distance = edges(ixedge) + tol/2;
        S.distance = S.distance(:);
        S.val      = lam;  
        if params.anisotropy            
            S.theta   = thetacents(ixtheta);
        end
    case 'cloud2'
        S.distance = iid(:,3);
        S.val      = lam;
        if params.anisotropy            
            S.theta   = thetacents(ixtheta);
        end
end

%% create plot if desired
if params.plotit
    switch params.type
        case {'default','gamma'}
            marker = 'o--';
        otherwise
            marker = '.';
    end
    
    if params.binomial
        marker='+--';
        'plot of binomial variogram : +'
    else
        marker=marker;
    end
    
    if ~params.anisotropy
        plot(S.distance,S.val,marker);
        axis([0 params.maxdist 0 max(S.val)*1.1]);
        xlabel('h');
        ylabel('\gamma (h)');
        title('(Semi-)Variogram');
    else
        [Xi,Yi] = pol2cart(repmat(S.theta,numel(S.distance),1),repmat(S.distance,1,numel(S.theta)));
        surf(Xi,Yi,S.val)
        xlabel('h y-direction')
        ylabel('h x-direction')
        zlabel('\gamma (h)')
        title('directional variogram')
%         set(gca,'DataAspectRatio',[1 1 1/30])
    end
end

