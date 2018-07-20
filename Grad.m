function [xGrad] = Grad(X)

    xGrad(:,:,1) = shiftAdj(X, [-1, 0], 'reflexive')-X;
    xGrad(:,:,2) = shiftAdj(X, [0, -1], 'reflexive')-X;
    
    
end