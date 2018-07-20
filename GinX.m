function [x] = GinX(x)
    a=max(max(x));
    b = min(min(x));
    if a> 255 || b<0
        x=10000000;
    else
        x=0;
    end
    
end