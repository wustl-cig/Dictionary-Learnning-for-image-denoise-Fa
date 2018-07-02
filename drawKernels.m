function drawKernels(d)
            %%% Visually illustrates dictionary patches
            %%%
            %%% Input:
            %%% - d: a set of 2D dictionary kernels
            %%%
            %%% U. S. Kamilov, CIG, WUSTL, 2017.
            
            [N,K] = size(d);
            kerSize = [sqrt(N),sqrt(N)];
            numKer = K;
            d = reshape(d,[kerSize,K]);
            d = padarray(d, [1, 1], min(d(:)), 'both');
            d = reshape(d, [kerSize+2, 1, numKer]);
            
            montage(d,...
                'Size', [NaN ceil(sqrt(numKer))],...
                'DisplayRange', [min(d(:)), max(d(:))]);
            colorbar;
 end