
function [u, warnmsg]=fem_eit_fwd_v2(numInc, backCond, cond, condOut, h, k, a, b, alpha, N, mypde)%, idx_circum, theta, mesh)
    %backCond, posNum, pos, xya, hex, hey, rad = phantom_infos;
    ccoefX = @(location, state)ccoeffunctionX(location, state, numInc, backCond, cond, h, k, a, b, alpha);
                                             %(location,state,backCond, posNum, pos, xya, hex, hey, rad)
    mygfun = @(location, state)mygfunN(location,state,N);   

    applyBoundaryCondition(mypde,"neumann", ...
                                 "Edge",1:mypde.Geometry.NumEdges, ...
                                 "g",mygfun);

    specifyCoefficients(mypde,"m",0,...
                              "d",0,...
                              "c",ccoefX,...
                              "a",0,...
                              "f",0);

    lastwarn('')%warning('');  % Clear last warning message
    results = solvepde(mypde);
    [warnmsg, msgid] = lastwarn;
    %sprintf(warnmsg)
    if not(isempty(warnmsg)) %then no need to do plots, so we exit function with return below
        u = 0;
        return
    end
    u = results.NodalSolution;
    

end
%%
function [X] = cart_ellipse(x, y, h, k, a, b, alpha)
        L = (((x-h).*cos(alpha) + (y-k).*sin(alpha)).^2)./(a^2);
        R = (((x-h).*sin(alpha) - (y-k).*cos(alpha)).^2)./(b^2);
        X = L+R;
end
%%
function gmatrix=mygfunN(location,state,N)    
    n1 = 1;
    nr = numel(location.x);
    gmatrix = zeros(n1,nr);
    r = sqrt(location.x^2 + location.y^2);
    if location.x > 0
        theta = atan(location.y/location.x);
    else 
        theta = atan(location.y/location.x) + pi;
    end   
    if N<=16
        gmatrix(1,:) = (1/sqrt(pi))*sin(N*theta);   
    else 
        gmatrix(1,:) = (1/sqrt(pi))*cos((N-16)*theta); 
    end
    
end
%%
function condOut=ccoeffunctionX(location, state, numInc, backCond, cond, h, k, a, b, alpha)
    condOut = ones(size(location.x))*backCond;
    x1 = location.x;
    x2 = location.y;
    for i = 1:numInc             
        X = cart_ellipse(x1, x2, h(i), k(i), a(i), b(i), alpha(i));

        condOut(X<=1)=cond(i);
    end
end