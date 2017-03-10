%%% functions needed to run variogramVOX_beta function

function iid = distmat2(X,dmax)

function iid = distmatsub1d(i) %au cas où il faut opitimiser distmat
    j  = (i+1:n)'; 
    d  = abs(X(i)-X(j));
    I  = d<=dmax;
    iid = [j(I) d(I)];
end

function iid = distmatsub2d(i)  %idem 2D
    j  = (i+1:n)'; 
    d = hypot(X(i,1) - X(j,1),X(i,2) - X(j,2));
    I  = d<=dmax;
    iid = [j(I) d(I)];
end
    
function iid = distmatsub(i) %idem 3D
    j  = (i+1:n)'; 
    d = sqrt(sum(bsxfun(@minus,X(i,:),X(j,:)).^2,2));
    I  = d<=dmax;
    iid = [j(I) d(I)];
end

% constrained distance function
%
% iid -> [rows, columns, distance]
n     = size(X,1);
nrdim = size(X,2);
if size(X,1) < 1000;
    [i,j] = find(triu(true(n)));%true:one matrix of dim n*n;triu:triagular upper part of the previous matrix
    %select unique pairs index (originally 2 sets of pair in a dist. matr
    if nrdim == 1;
        d = abs(X(i)-X(j));
    elseif nrdim == 2;
        d = hypot(X(i,1)-X(j,1),X(i,2)-X(j,2));
    else
        d = sqrt(sum((X(i,:)-X(j,:)).^2,2));%3D case distance computation
    end
    I = d<=dmax;
    iid = [i(I) j(I) d(I)];%table with pairs index ans associated distance by row
else
    ix = (1:n)';%idem optimized
    if nrdim == 1;
        iid = arrayfun(@distmatsub1d,(1:n)','UniformOutput',false);
    elseif nrdim == 2;
        % if needed change distmatsub to distmatsub2d which is numerically
        % better but slower
        iid = arrayfun(@distmatsub,(1:n)','UniformOutput',false);
    else
        iid = arrayfun(@distmatsub,(1:n)','UniformOutput',false);
    end
    nn  = cellfun(@(x) size(x,1),iid,'UniformOutput',true);  
    I   = nn>0;
    ix  = ix(I);
    nn  = nn(I);
    nncum = cumsum(nn);
    c     = zeros(nncum(end),1);
    c([1;nncum(1:end-1)+1]) = 1;
    i = ix(cumsum(c));
    iid = [i cell2mat(iid)];
    
end
end