function [dsig,hist] = recon_sparsity(invpm,fem,verbose,sigt)

dsig = zeros(fem.sdof,1); % initial inhomog.
hist   = [];                % iter. history
MaxF   = 5;

[gradfcn,objfcn,regval] = fcnl_eval(invpm,fem,dsig);
for iter = 1:invpm.maxit
    % Sobolev gradient
    ffs = sctermpw(fem,gradfcn);
    [gradfcn] = dirichsolve(invpm.KKs,zeros(length(fem.bcdof),1),ffs,fem);
    % Barzilai-Borwein step length with safeguard
    if iter==1
        dscdir = -gradfcn;
        dscdir0 = -gradfcn;
        step_bb = 0;
        step = 1;
        dsig0 = dsig;
        dsig = dsig + step*dscdir;
        dsig = sign(dsig).*max(abs(dsig)-step*invpm.eta,0);
        dsig = max(dsig,-.99);
        [gradfcn,objfcn,regval] = fcnl_eval(invpm,fem,dsig);
    else
        dscdir = -gradfcn;
        step_bb = -((dsig-dsig0)'*invpm.KKs*(dscdir-dscdir0))/((dscdir-dscdir0)'*invpm.KKs*(dscdir-dscdir0));
        step = min(max(step_bb,1),2048);
        dscdir0 = dscdir;
    end

    % implement safe-guard rule
    if step_bb>0 % dscdir: descent direction
        dsig0 = dsig;
        [step,gradfcn,objfcn,regval,dsig] = bb_sgd(invpm,fem,dsig,fcnvalmax,dscdir,step);
    elseif ((step_bb<0) & iter >1)
        dsig0= dsig;
        dsig = dsig + step*dscdir;
        dsig = sign(dsig).*max(abs(dsig)-step*invpm.eta,0);
        dsig = max(dsig,-.99);
        [gradfcn,objfcn,regval] = fcnl_eval(invpm,fem,dsig);
    end

    if mod(iter,100)==0
        filename  = '';
    end

    % record & display results
    if exist('sigt')
        err = [0];%%%%%(dsig+invpm.sig0-sigt)'*invpm.MM*(dsig+invpm.sig0-sigt);
    else
        err = [0];
    end
    hist = [hist;iter objfcn regval objfcn+regval step_bb step err];
    
    ind = max(1,iter-MaxF+1):iter;
    fcnvalmax = max(hist(ind,4));

    if verbose
        figure(3), clf
        subplot(2,1,1), pdeplot(fem.p,fem.e,fem.t,'xydata',invpm.sig0+dsig), drawnow, colormap('hsv')
        subplot(2,1,2), pdecont(fem.p,fem.t,dsig,20), drawnow, colorbar
        display(sprintf('iter=%0.5g, objval=%0.5g, regval=%0.5g, fcnt=%0.5g, step=%0.5g, step_bb = %0.5g, fcnvalmax = %0.5g', ...
            iter,objfcn,regval,objfcn+regval,step,step_bb,fcnvalmax))
    end
    
    % stopping criterion to be implemented
    if abs(step) < 1e-3
        break
    end
end