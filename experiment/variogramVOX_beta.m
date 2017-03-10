function S = variogramVOX_beta(x,y,binomial_var,varargin)

% isotropic experimental ordinary or binomial(semi-)variogram

%%% use distmatBound.m and parseargs.m to run

% Syntax:
%   d = variogramVOX(x,y,binomial_var)
%   d = variogramVOX(x,y,binomial_var'propertyname','propertyvalue',...)
%   NB: in ordinary variogram, binomial_var is not used but has to be set
%   to a non-empty value, e.g. 0

%Input:
%   x - array with coordinates. Each row is a location in a 
%       size(x,2)-dimensional space (e.g. [x y elevation])
%
%   y - column vector with values of the locations in x. 
%
%   binomial_var - column vector with number of experiment for this row
%
%
% Propertyname/-value pairs:
%   nrbins - number bins the distance should be grouped into
%            (default = 20)
%   maxdist - maximum distance for variogram calculation
%            (default = maximum distance in the dataset / 2)
%   type -   'gamma' returns the variogram value (default)
%            'cloud1' returns the binned variogram cloud
%            'cloud2' returns the variogram cloud
%   plotit - true -> plot variogram
%            false -> don't plot (default)
%   subsample - number of randomly drawn points if large datasets are used.
%               scalar (positive integer, e.g. 3000)
%               inf (default) = no subsampling
%   INACTIVE anisotropy - false (default), true (works only in two dimensions)
%   thetastep - if anisotropy is set to true, specifying thetastep 
%            allows you the angle width (default 30°)
%   
%   
% Output:
%   d - structure array with distance and gamma - vector

%10/03/2017
% Max, from Wolfgang Schwanghart

%Experimental variogram computing

% % Exemple Arguments settings for test
% x=[RDI(:,1),RDI(:,2),RDI(:,3)];%matrice de coordonnées des points observés
% y=RDI(:,4);%variable RDI
% binomial_var=RDI(:,5);%variable N=Nt-Nb
% 
% nrbins=15;
% maxdist=0.35;
% plotit=true;
% binomial=1;

%varargin={'nrbins',nrbins,'maxdist',maxdist,'binomial',binomial,'plotit',plotit};%list of non-default arguments values

% extent of dataset
minx = min(x,[],1);
maxx = max(x,[],1);
maxd = sqrt(sum((maxx-minx).^2));
nrdims = size(x,2);

% check input using PARSEARGS (parameters values)
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



% size checking
if size(y,1) ~= size(x,1);
    error('x and y must have the same number of rows')
end
if params.binomial
   if size(binomial_var,1) ~= size(y,1);
    error('y and binomial variable must have the same number of rows')
   end
end 

% check for nans
if params.binomial
    II = any(isnan(x),2) | isnan(y) | isnan(binomial_var);%tests matrixs rows
    x(II,:) = [];
    y(II)   = [];
    binomial_var(II) = [];
    if unique(II) ~= 0;
        warning('NaN values have been removed')
    end
else
    II = any(isnan(x),2) | isnan(y);%tests matrixs rows
    x(II,:) = [];
    y(II)   = [];
    if unique(II) ~= 0;
        warning('NaN values have been removed')
    end
end

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





% if specified, take only a subset of the data at random;
if ~isinf(params.subsample) && numel(y)>params.subsample;
    IX = randperm(numel(y),params.subsample);
    x  = x(IX,:);
    y  = y(IX,:);
end

% calculate bin tolerance
tol = params.maxdist/params.nrbins; %pas de taille des classes

% calculate constrained distance matrix
iid = distmatBound(x,params.maxdist);

% calculate squared difference between values of coordinate pairs BINOMIAL MODIFIED
if params.binomial
    Npair_product=binomial_var(iid(:,1)).*binomial_var(iid(:,2));
    Npair_sum=binomial_var(iid(:,1))+binomial_var(iid(:,2));
    ind=find(binomial_var>=30);%TO INCREASE TO BE MORE CONSERVATIVE
    expct=mean(y(ind));%OBSERVED expectation of RDI
    var2=var(y(ind));%OBSERVED variance of RDI
    Nterm=Npair_product./Npair_sum;
    
        lam = (Nterm.*((y(iid(:,1))-y(iid(:,2))).^2)-expct*(1-expct)+var2);%/!\ valeurs légèrement négative
else
        lam      = (y(iid(:,1))-y(iid(:,2))).^2;
end
    iid=[iid lam];

% calculate variogram
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

% create plot if desired
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
end