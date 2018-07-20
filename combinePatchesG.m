function x = combinePatchesG(g,imgSize,DC)
            
            [N,P,d] = size(g);
            kerSize = [sqrt(N),sqrt(N)];
            osize=[imgSize,d];
             x = zeros(osize);
            idMax = imgSize-kerSize+1;
            g = g+ DC;
            for k = 1:2
                count = 1;
                for i = 1:idMax(1)
                    for j = 1:idMax(2)
                        x(i:i+kerSize(1)-1,j:j+kerSize(2)-1,k) = x(i:i+kerSize(1)-1,j:j+kerSize(2)-1,k)+ reshape(g(:,count,k),kerSize) ;
                        count = count +1;
                    end
                end
                x = 1/N * x;
            end
end