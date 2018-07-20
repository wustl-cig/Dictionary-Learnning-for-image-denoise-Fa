function [x] = Gtran(Grad)
            x1 = Grad(:,:,1);
            x2 = Grad(:,:,2);
            x1 = shiftAdj(x1, [-1,0], 'reflexive')-x1;
            x2 = shiftAdj(x2, [0,-1], 'reflexive')-x2;
            
            x = x1 + x2;