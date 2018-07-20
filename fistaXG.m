function [X] = fistaXG(X,Y,A,alpha,lam)


fprintf("fistaX");
% atran = @(z) A'*z;
[N,k]=size(A);
[a,b]=size(Y);
[k,p]=size(alpha);

% useful shortcuts


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% optimize: gradient method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    xhat0 = X; % initial value
%     gamma = 1; % step-size
    numIter = 20;
    beta = 0.5;
    
    backTol = 1e-20;
    % initialize
    xhat = xhat0;
    xhatprev = xhat;
    tprev = 1;

    % differentiable function

    f = @(x) sum(vecnorm(x-Y).^2) + lam*sum(sum(vecnorm(ExpatchG(sqrt(N),Grad(x)) - reshape(A*alpha,[N,p/2,2])).^2));
    grad_f = @(x) 2*(x-Y) + 2*lam*N*Gtran(combinePatchesG(ExpatchG(sqrt(N),Grad(x))-reshape(A*alpha,[N,p/2,2]),[a,b],0));


    % proximable function
    g = @(x) GinX(x);
    prox_g = @(x) 0.*(x<0)+255.*(x>255)+x.*(x>=0&x<=255);

    % complete function
    h = @(x) f(x) + g(x);
    H = zeros(1,numIter);

%     figure(102)
        for iter = 1:numIter
            
            gamma = 1; %initial step-size
            % update theta
            t = 0.5*(1+sqrt(1+4*tprev*tprev));
    
            % acceleration
            s = xhat + ((tprev-1)/t)*(xhat-xhatprev);
    
         % compute the next iterate
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
            % begin print
            fprintf("fista X with step" + gamma+",");
             temp1 = f(xhat);
            fprintf("[iteration"+ iter+"/" +numIter+"]" );
            fprintf("f:"+ f(xhat)+".");
            fprintf("g:"+ g(xhat)+".");
            fprintf("f+g:"+ h(xhat)+".\n");
        
             H(1,iter) = h(xhat);
        end
       
%        plot(1:iter,H);
%        drawnow;



X=xhat;
end