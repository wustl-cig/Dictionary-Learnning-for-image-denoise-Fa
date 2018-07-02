function [alpha] = fista(Old,X,A,rho,tau)

% alpha is the sparse coding of X based on dictionary A and premeter lam.
fprintf("fista");
amult = @(x) A*x;
atran = @(z) A'*z;
[N,k]=size(A);
[N,p]=size(X);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% optimize: gradient method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    xhat0 = zeros(k, p); % initial value
    gamma = 0.1; % step-size
    numIter = 40;
    beta = 0.5;
    
    backTol = 1e-12;
    % initialize
    xhat = xhat0;
    xhatprev = xhat;
    tprev = 1;

    % differentiable function
    % f = @(x) 0.5*norm(X-amult(x), 'fro')^2;
    f = @(x) rho*sum(vecnorm(X-amult(x)).^2);
    % grad_f = @(x) sum(atran(amult(x)-X),2);
    grad_f = @(x) rho*2*atran(amult(x)-X);

    % proximable function
    g = @(x) tau*norm(x(:),1);
     prox_g = @(x) wthresh(x,'s',tau*gamma);
%     prox_g = @(x) getAlpha(x,tau/rho,gamma);
    % complete function
    h = @(x) f(x) + g(x);
    H = zeros(1,numIter);



        for iter = 1:numIter
            gamma = 0.1; %initial step-size
           
            
            % update theta
            t = 0.5*(1+sqrt(1+4*tprev*tprev));
    
            % acceleration
            s = xhat + ((tprev-1)/t)*(xhat-xhatprev);
            
             % compute the next iterate
            xhatnext = prox_g(s - gamma*grad_f(s));
            
            % backtracking line search
            tol = f(s)-f(xhatnext)...
                +(0.5/gamma)*(sum(vecnorm(xhatnext-s+gamma*grad_f(s)).^2))...
                -(0.5*gamma)*sum(vecnorm(grad_f(s)).^2);
            
            while (tol<0 && backTol>0)
%                 fprintf("shrink!")
                gamma = gamma*beta;
               
%                 t = theta*gamma*sum(vecnorm(grad_f(s)).^2);
                tol = f(s)-f(xhatnext)...
                    +(0.5/gamma)*(sum(vecnorm(xhatnext-s+gamma*grad_f(s)).^2))...
                    -(0.5*gamma)*sum(vecnorm(grad_f(s)).^2);                  % recompute the next iterate
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
            fprintf("fista alpha with step" + gamma+",");
            fprintf("[iteration"+ iter+"/" +numIter+"]" );
            fprintf("f:"+ f(xhat)+".");
            fprintf("g:"+ g(xhat)+".");
            fprintf("f+g:"+ h(xhat)+".\n");
            H(iter) = h(xhat);
            
        end
       



alpha=xhat;
end
