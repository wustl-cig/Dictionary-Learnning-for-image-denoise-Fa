function x = combinePatches(g,imgSize,DC)
            
            [N,P] = size(g);
            kerSize = [sqrt(N),sqrt(N)];
            x = zeros(imgSize);
            idMax = imgSize-kerSize+1;
            g = g+ DC;
            count = 1;
            for i = 1:idMax(1)
                for j = 1:idMax(2)
                    x(i:i+kerSize(1)-1,j:j+kerSize(2)-1) = x(i:i+kerSize(1)-1,j:j+kerSize(2)-1)+ reshape(g(:,count),kerSize) ;
                    count = count +1;
                end
            end
            x = 1/N * x;
        end