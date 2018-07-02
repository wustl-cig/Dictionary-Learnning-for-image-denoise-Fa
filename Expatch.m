function [X]=Expatch(n,image)
       kerSize = [n,n];
       imgSize = size(image);
       idMax= imgSize-kerSize+1;
       X = zeros(n*n,idMax(1,1)*idMax(1,2));
       count = 1;
       for i = 1:idMax(1)
           for j = 1:idMax(2)
               gTemp = image(i:i+kerSize(1)-1,j:j+kerSize(2)-1);% row iters
               X(:,count)= reshape(gTemp,[],1);
               count = count+1;
           end
       end
            % whitening
%        DC = mean(X,1);
%        X = X-DC; 
       
        
end