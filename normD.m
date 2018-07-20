function [D] = normD(D)
    [N,k]= size(D);
    for o=1:k 
        if norm(D(:,o))>1
           
             D(:,o)=D(:,o)/norm(D(:,o));
        end
  
    end 
%     D./repmat(sqrt(sum(x.^2+0.0001,1)),size(x,1),1);

end