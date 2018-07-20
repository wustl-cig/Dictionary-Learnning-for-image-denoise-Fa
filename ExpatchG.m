function [GradX] = ExpatchG(n,I)
     GradX(:,:,1)=Expatch(n,I(:,:,1));
     GradX(:,:,2)=Expatch(n,I(:,:,2));

end