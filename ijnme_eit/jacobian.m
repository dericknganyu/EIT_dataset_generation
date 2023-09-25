function [J,FN0d] = jacobian(fem,sigma0,flux)
%------------------------------------------------------------------------%
% Purpose: computing linearized inverse problem for electrical impedance %
% tomography by considering NtD map, sensitivity analysis                %
% written by Bangti Jin (btjin@math.uni-bremen.de) May 5, 2010           %
% Input: p,e,t:  finite-element parameter                                %
%        sigma0: background conductivity                                 %
%        sigmat: conductivity with inclusions                            %
%        flux:   input current                                           %
%------------------------------------------------------------------------%

bcdof = fem.bcdof;
Ned   = length(bcdof);

%----------- stiffness matrix & source term ------------------------%
KK0 = stiffmat(fem,sigma0);
[N] = size(flux,2);
FN = [ ];
FNn = [ ];
FN0d =[ ];
for i = 1:N
    % approximation
    [u,ud,un] = neumansolve(KK0,flux(:,i),zeros(fem.sdof,1),fem);
    FN   = [FN u];
    FN0d = [FN0d;ud];
    FNn  = [FNn un];
end

J = [ ];    % change only in interior part of domain
for j=1:N
    vec = [];
    for inode = 1:fem.nnode
        if isempty(find(bcdof==inode))  % not on boundary
            ff        = fstermspars(fem,inode,fem.t,FN(:,j));
            [u,ud,un] = neumansolve(KK0,zeros(Ned,1),ff,fem);
            vec       = [vec ud];
        end
    end
    display(['flux computed: ' num2str(j)])
    J = [J; vec];
end