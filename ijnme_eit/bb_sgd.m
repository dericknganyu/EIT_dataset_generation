function [step,gradfcn,objfcn,regval,dsig] = bb_sgd(invpm,fem,dsig,fcnvalmax,dscdir,sinit)
%--------------------------------------------------------------------%
% implementing safeguard strategy for Barzilai-Borwein rule          %
% Created on: Sept. 18, 2010 by  Bangti Jin (btjin@math.tamu.edu)    %
%--------------------------------------------------------------------%

epsilon = 1e-5;
etaS    = 2;
eta     = invpm.eta;

step = sinit*etaS;
for it = 1:100
    step = step / etaS;
    % form gradient of data-fitting function
    dsigtp = dsig + step*dscdir;
    dsigtp = sign(dsigtp).*max(abs(dsigtp)-step*eta,0);
    dsigtp = max(dsigtp,-.99);
    iterdiff = epsilon*step/2* (dsigtp - dsig)'*invpm.KKs*(dsigtp-dsig);
    [gradfcn,objfcn,regval] = fcnl_eval(invpm,fem,dsigtp);
    fcnval = objfcn + regval;
    if (fcnval <= fcnvalmax - iterdiff)|(abs(step)<=1e-3)
        break
    end
end
dsig = dsigtp;