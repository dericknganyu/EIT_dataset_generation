function [gradfcn,objfcn,regval] = fcnl_eval(invpm,fem,dsig)
%-------------------------------------------------------------------------%
% Purpose: calculate gradient, objective functional value, reg. term etc. %
% Copyright: Bangti Jin (kimbtsing@yahoo.com.cn)                          %
% Created at Sept. 18, 2010                                               %
%-------------------------------------------------------------------------%

% extracting parameters
Neig = invpm.Neig;
sig0 = invpm.sig0;
flux = invpm.flux;
temp = invpm.temp;
wt   = invpm.wt;
eta  = invpm.eta;

% Neig,MM,sigma0,dsig,p,e,t,flux,temp,wt,eta)
KKi     = stiffmat(fem,sig0+dsig);
ffo     = zeros(fem.sdof,1);
gradfcn = zeros(1,fem.nel);
objfcn  = 0;
for j = 1 : Neig
    % solving forward problem
    [uN,uNt,uNn] = neumansolve(KKi,flux(:,j),ffo,fem);
    rd  = uNt - temp(:,j);
    % solving adjoint problem
    [uA,uAt,uAn] = neumansolve(KKi,rd,    ffo,fem);
    % objective function value
    objfcn = objfcn + wt(j) * rd'*invpm.MMb*rd/2;
    % form gradient of data-fitting function
    [uNx,uNy] = pdegrad(fem.p,fem.t,uN);
    [uAx,uAy] = pdegrad(fem.p,fem.t,uA);
    gradfcn =  gradfcn - wt(j) * (uNx.*uAx + uNy.*uAy);
end
regval = sum(invpm.MM*abs(dsig))*eta;