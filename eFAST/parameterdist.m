%% Parameter distributions for sampling scheme 
% [pmin pmax]: min and max values of the range of variation
% [pmean pstd]: mean and standard deviation for distributions other than uniform (normal and lognormal)
% type: type of distributions (pdf). 
% Uniform ['unif'], nomral ['norm'], triangular

function Xdist = parameterdist(X,pmax,pmin,pmean,pstd,type)
Xdist = NaN(size(X));
[~, nparam] = size(Xdist);
switch lower(type)
    case {'unif'}
        for k=1:nparam %loop through parameters
            Xdist(:,k)=(X(:,k).*(pmax(k)-pmin(k)))+pmin(k);
        end
    case {'triangular'}
        for k=1:nparam %loop through parameters
            aa = min(pmin(k),pmax(k));
            cc = max(pmin(k),pmax(k)); 
            pd = makedist('Triangular','A',aa,'B',pmean(k),'C',cc);  % bound [A,C], peak at B
            Xdist(:,k) = icdf(pd,X(:,k));
        end
    case {'norm'}
        for k=1:nparam %loop through parameters
            Xdist(:,k) = norminv(X(:,k),pmean(k),pstd(k));
        end
%     case {'lognorm'}
%         for k=1:nparam %loop through parameters
%             y = logninv(X(:,k),pmean(k),pstd(k));
%             Xdist(:,k) = exp(y);
%             keyboard
%         end
    otherwise
        disp('Unknown pdf')
end