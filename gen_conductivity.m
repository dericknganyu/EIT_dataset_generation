function [numInc, backCond, cond, condOut, h, k, a, b, alpha, kx, ky] = gen_conductivity(mesh, max_numInc, texture)
    nodes = mesh.Nodes;
    
    x1 = nodes(1,:);
    x2 = nodes(2,:);
    
    numInc = randi(max_numInc);

    cond = zeros(1, numInc);
    kx   = zeros(1, numInc);
    ky   = zeros(1, numInc);
    
    
    backCond=1;
    condOut = ones(size(x1))*backCond;
    
    [h, k, a, b, alpha] = sampleInclusions(numInc);
    
    if islogical(texture)
        if texture
            for i = 1:numInc             
                X = cart_ellipse(x1, x2, h(i), k(i), a(i), b(i), alpha(i));
                X1 = x1(X<=1);  X2 = x2(X<=1);
                kx(i) = unifrnd(5, 15, 1, 1);%unifrnd(1, 2, 1, 1); % 
                ky(i) = unifrnd(5, 15, 1, 1);%unifrnd(1, 2, 1, 1); %
                centre = [h(i), k(i)];
                res = 0.5*(1 + add_texture(X2, X1, kx(i), ky(i), alpha(i), centre)); %scaling so that [0, 1]
                
                cond_opt  = [0, 1]; 
                cond_idx  = randi([1, 2], 1); %chooseing between 0 and 1
                cond(i) =  cond_opt(cond_idx);
                
                res = 0.6*res + 0.2 + cond(i)*(0.2*res + 1); %if cond(i) was 0 it reduces to 0.6*res + 0.2 in [0.2, 0.8]
                                                            %if cond(i) was 1 it reduces to 0.8*res + 1.2 in [0.8, 2.0]  
                
                condOut(X<=1)=res;% cond(i);
            end
        else
            for i = 1:numInc             
                cond_opt  = [rand*0.6+0.2, rand*0.8 + 1.2];  %sigma in [0.2,0.8] or [1.2, 2]
                cond_idx  = randi([1, 2], 1);
                cond(i) =  cond_opt(cond_idx);
                X = cart_ellipse(x1, x2, h(i), k(i), a(i), b(i), alpha(i));

                condOut(X<=1)=cond(i);
            end
        end       
    else
        % Do nothing and we have the case of constant conductivity in all domain        
        fprintf('\nTexture was neither true or false')
        fprintf('\ndefaulting to case of constant conductivity of 1\n\n')
    end 
end
%%
function [z] = add_texture(x, y, kx, ky, angle, centre)

    x_rot = centre(1) + x*cos(angle) - y*sin(angle);
    y_rot = centre(2) + x*sin(angle) + y*cos(angle);

    z = 0.5*(sin(kx*x_rot) + sin(ky*y_rot)); %0.5*(2 + sin(kx*x_rot) + sin(ky*y_rot)); %[-1, 1] #Delete that 2!!
    
end
%%
function [X] = cart_ellipse(x, y, h, k, a, b, alpha)
        L = (((x-h).*cos(alpha) + (y-k).*sin(alpha)).^2)./(a^2);
        R = (((x-h).*sin(alpha) - (y-k).*cos(alpha)).^2)./(b^2);
        X = L+R;
end
%% 
function [h, k, a, b, alpha] = sampleInclusions(numInc, test_step, tolerance)
    narginchk(1, 3);
    if nargin < 3
        tolerance = 5;%200;
    end
    if nargin < 2
        test_step = 200;
    end
    
    h = zeros(1, numInc); k = zeros(1, numInc);
    a = zeros(1, numInc); b = zeros(1, numInc);
    alpha = zeros(1, numInc); 
  

    [h(1), k(1), a(1), b(1), alpha(1)] = sampleEllipse(test_step);

    %tol = 0;    
    for i = 2:numInc 
        overlap = true;         
        tol = 0;
        while overlap && tol<tolerance
            overlap = false;
            [h(i), k(i), a(i), b(i), alpha(i)] = sampleEllipse(test_step);
            for j = 1:i-1  
                if norm([h(j), k(j)] - [h(i), k(i)]) < a(j) + a(i)
                    overlap = true; % means generated inclusion 'overlaps' with another
                    break             % exit for loop and generate new sample
                end %if norm([h(j), k(j)] - [h(i), k(i)]) < a(j) + a(i); overlap = true; end
            end
            tol = tol + 1; %when tol finally equals tolerance, the while loop ends
                           %Meaning we have tried to sampling non-overlapping
                           %ellipses more than we can tolerate and they keep overlapping, so we should restart
                           %sampling of the 1st sample
        end
        if tol == tolerance%; break; end     
            [h, k, a, b, alpha] = sampleInclusions(numInc, test_step, tolerance);
            break;
        end
        
    end
     
end

function [h, k, a, b, alpha] = sampleEllipse(test_step, tolerance)
    %# (h, k) centre of ellipse
    %# (a, b) major and minor axis of ellipse
    %# alpha angle of rotation
    %# test_step: just a variable needed to identify the points on an ellpise via their polar coordinate so we can check if it is within the unit radius
    %# tolerance: how many times we are allowed to sample b wrongly before restarting the sampling process of h, k, alpha, a

    narginchk(1, 2);
    if nargin<2,  tolerance=50; end
    
    h = rand*(1.8) - 0.9; %[-0.9, 0.9]
    k = rand*(1.8) - 0.9; %[-0.9, 0.9]
    alpha = rand*2*pi; %[0, 2pi] 
    a = rand*(0.8) + 0.1;   %[0.1, 0.9] 
    b = rand*(a-0.1) + 0.1; %[0.1, a]   

    theta = linspace(0,2*pi,test_step);
    x = h + a*cos(alpha).*cos(theta) - b*sin(alpha).*sin(theta);
    y = k + a*sin(alpha).*cos(theta) + b*cos(alpha).*sin(theta);
    i = 0;
    while any(x.^2 + y.^2 > 0.9) 
        b = rand*(a-0.1) + 0.1; %[0.1, a] %[0.3, 0.5]
        x = h + a*cos(alpha).*cos(theta) - b*sin(alpha).*sin(theta);
        y = k + a*sin(alpha).*cos(theta) + b*cos(alpha).*sin(theta);
        if i == tolerance
            i = 0;
            [h, k, a, b, alpha] = sampleEllipse(test_step);
        end
        i = i+1;
    end
end

%%
% function plot_Elipse_in_domain(x, y)
%     %plotting circle
%     th = 0:pi/100:2*pi;
%     xunit = cos(th);
%     yunit = sin(th);
%     hold on
%     plot(xunit, yunit)
%     numInc = size(x, 1);
%     for i=1:numInc
%         plot(x(i, :), y(i, :),'.-')
%     end
%     axis equal
%     axis off
%     hold off
% end
%%
%%
% function [x, y] = calculateEllipse(h, k, a, b, alpha, steps)
%     %# This functions returns points to draw an ellipse
%     %#
%     %#  @param h     X coordinate
%     %#  @param k     Y coordinate
%     %#  @param a     Semimajor axis
%     %#  @param b     Semiminor axis
%     %#  @param angle Angle of the ellipse (in degrees)
%     %#
%     narginchk(5, 6);
%     if nargin<6, steps = 100; end

%     theta = linspace(0, 2*pi, steps);%' .* (pi / 180);

%     x = h + a*cos(alpha).*cos(theta) - b*sin(alpha).*sin(theta);
%     y = k + a*sin(alpha).*cos(theta) + b*cos(alpha).*sin(theta);

%     if nargout==1, [x, y]; end
% end
%%

