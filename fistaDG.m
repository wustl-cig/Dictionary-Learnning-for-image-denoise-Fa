function [D] = fistaDG(X,alpha,Dictionary,rho)
fprintf("fistaD")
 [k,p,d]=size(alpha);
 [N,p]=size(X);
% [a,b]=size(Y);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% minimize f(D)=sum_i=1^p ||R_i*X-D*alpha_i|| + 11_D(c)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% atran = @(z) A'*z;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% optimize: gradient method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    xhat0 = Dictionary; % initial value
%     gamma = 0.001; % step-size
    numIter = 50;
    beta = 0.5;
    backTol = 1e-12;
    % initialize
    xhat = xhat0;
    xhatprev = xhat;
    tprev = 1;

    % differentiable function
   
     f = @(x) rho * sum(vecnorm(X-x*alpha).^2);
     grad_f=@(x) 2*rho*(x*alpha-X)*alpha';

    

    % proximable function
    g = @(x) 0;
    prox_g = @(x) normD(x);

    % complete function
    h = @(x) f(x) + g(x);
    H = zeros(1,numIter);


        for iter = 1:numIter
            gamma = 1; %initial step-size
            % update theta
            t = 0.5*(1+sqrt(1+4*tprev*tprev));
    
            % acceleration
            s = xhat + ((tprev-1)/t)*(xhat-xhatprev);
    
         % compute the next iterate
%             temp1 = gamma*grad_f(s);
            xhatnext = prox_g(s - gamma*grad_f(s));
            
            % backtracking line search
            tol = f(s)-f(xhatnext)+0.5/gamma*(sum(vecnorm(xhatnext-s+gamma*grad_f(s)).^2))-0.5*gamma*sum(vecnorm(grad_f(s)).^2);
           
            while tol<0
%                 fprintf("shrink!")
                gamma = gamma*beta;
               
%                 t = theta*gamma*sum(vecnorm(grad_f(s)).^2);
                tol = f(s)-f(xhatnext)+0.5/gamma*(sum(vecnorm(xhatnext-s+gamma*grad_f(s)).^2))-0.5*gamma*sum(vecnorm(grad_f(s)).^2);                  % recompute the next iterate
                xhatnext = prox_g(s - gamma*grad_f(s));
                 if gamma < backTol
                    gamma = backTol;
                    break;
                end
            end
    
            % update the iterate
            xhatprev = xhat;
            xhat = xhatnext;
            tprev = t;
            
            
            H(1,iter) = h(xhat);
            fprintf("fista Dictionary with step" + gamma+",");
%             fprintf("fista Dictionary")
            fprintf("[iteration"+ iter+"/" +numIter+"]" );           
            fprintf("f+g:"+ h(xhat)+".\n");
            
        end
%   plot(H);
%   drawnow;



D=xhat;
end