function [Q]=orthtrans(F,type,KK,fem)
%-------------------------------------------------------------------------%
% perform Gram-Schmidt process for computing orthonormal basis of vectors %
% F with the inner product given in type                                  %
% JIN Bangti (kimbtsing@yahoo.com.cn) Feb. 24, 2009                       %
%-------------------------------------------------------------------------%

gcoord   = fem.gcoord;
sdof     = fem.sdof;
[n,Neig] = size(F);
Q        = [];

switch type
    case 'L2'
        % compute the matrix for inner product first
        MM = massmatbdy(fem); 
        for i=1:Neig
            q=F(:,i);
            for j=1:i-1
                q = q - (q'*MM*Q(:,j))*Q(:,j);
            end
            q = q/sqrt(q'*MM*q);
            Q = [Q q];
        end
    
    case 'Hm'
        % compute the inner product on the fly
        uN = [];
        for i=1:Neig
            q = F(:,i);
            for j=1:i-1
                [uNq] = neumansolve(KK,q,zeros(sdof,1),fem);
                    q = q - (uNq'*KK*uN(:,j))*Q(:,j);
            end
            [uNq] = neumansolve(KK,q,zeros(sdof,1),fem);
            cnst  = sqrt(uNq'*KK*uNq); % normalizing constant
            q     = q/cnst;            % normalization
            uN    = [uN uNq/cnst];     % normalized solution for saving computing in inner loop
            Q     = [Q q];
        end
    otherwise
        error('Undefined inner product!')
        
end