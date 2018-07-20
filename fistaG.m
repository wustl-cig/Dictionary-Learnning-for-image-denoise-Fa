function [alpha] = fistaG(Old,X,A,rho,tau)

% alpha is the sparse coding of X based on dictionary A and premeter lam.
fprintf("fista alpha");
[N,k,d]=size(A);
[N,p,d]=size(X);
amult = @(x) reshape([A(:,:,1)*x(:,:,1),A(:,:,2)*x(:,:,2)],[N,p,2]);
atran = @(z) [A(:,:,1)'*z(:,:,1),A(:,:,2)'*z(:,:,2)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% optimize: gradient method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    xhat0 = zeros(k, p,2); % initial value
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
    f = @(x) rho*sum(vecnorm([X(:,:,1),X(:,:,2)]-amult(x)).^2);
    % grad_f = @(x) sum(atran(amult(x)-X),2);
 
    grad_f = @(x) rho*2*atran(amult(x)-[X(:,:,1),X(:,:,2)]);

    % proximable function
    g = @(x) tau*norm(x(:),1);
     prox_g = @(x) wthresh(x,'s',tau*gamma);
%     prox_g = @(x) getAlpha(x,tau/rho,gamma);
    % complete function
    h = @(x) f(x) + g(x);
    H = zeros(1,numIter);



        for iter = 1:numIter
            gamma = 0.1; %initial step-size
           
            
            % data processing
            
            
            % update theta
            t = 0.5*(1+sqrt(1+4*tprev*tprev));
    
            % acceleration
            s = xhat + ((tprev-1)/t)*(xhat-xhatprev);
            
             % compute the next iterate
            xhatnext = prox_g([s(:,:,1),s(:,:,2)] - gamma*grad_f(s));
            
            % backtracking line search
            tol = f([s(:,:,1),s(:,:,2)])-f(xhatnext)...
                +(0.5/gamma)*(sum(vecnorm(xhatnext-[s(:,:,1),s(:,:,2)]+gamma*grad_f(s)).^2))...
                -(0.5*gamma)*sum(vecnorm(grad_f(s)).^2);
            
            while (tol<0 && backTol>0)
%                 fprintf("shrink!")
                gamma = gamma*beta;
               
%                 t = theta*gamma*sum(vecnorm(grad_f(s)).^2);
               tol = f([s(:,:,1),s(:,:,2)])-f(xhatnext)...
                +(0.5/gamma)*(sum(vecnorm(xhatnext-[s(:,:,1),s(:,:,2)]+gamma*grad_f(s)).^2))...
                -(0.5*gamma)*sum(vecnorm(grad_f(s)).^2);
                xhatnext = prox_g([s(:,:,1),s(:,:,2)] - gamma*grad_f(s));
                 if gamma < backTol
                    gamma = backTol;
                    break;
                end
            end
    
            
            xhatnext=reshape(xhatnext,[k,p,2]);
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
