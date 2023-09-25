function [cnst]=orthtransnorm(F,type,KK,fem)
%-------------------------------------------------------------------------%
% perform Gram-Schmidt process for computing orthonormal basis of vectors %
% F with the inner product given in type                                  %
% JIN Bangti (kimbtsing@yahoo.com.cn) Feb. 24, 2009                       %
%-------------------------------------------------------------------------%

gcoord   = fem.gcoord;
sdof     = fem.sdof;
[n,Neig] = size(F);

cnst = [];
switch type
    case 'L2'
        % compute the matrix for inner product first
        MM = massmatbdy(fem); 
        for i=1:Neig
            q=F(:,i);
            cnst = [cnst sqrt(q'*MM*q)];
        end
    
    case 'Hm'
        % compute the inner product on the fly
        uN = [];
        for i=1:Neig
            q = F(:,i);
            [uNq] = neumansolve(KK,q,zeros(sdof,1),fem);
            cnst = [cnst sqrt(uNq'*KK*uNq)]; % normalizing constant
        end
        
    otherwise
        error('Undefined inner product!')        
end