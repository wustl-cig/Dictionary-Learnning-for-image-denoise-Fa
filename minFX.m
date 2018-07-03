clear all;

 Io=double(imread('data/barbara.jpg'))/255;
 
sigma=20/255;
[a,b]=size(Io);
n=9;
 I=Io+(sigma)*randn(size(Io));
X=Expatch(n,I);
me=mean(X);
[N,p]=size(X);
K=256;
Dictionary=zeros(n*n,K); 

for k=0:1:K-1 
    test=X(:,randi(p));
%    	V=dct2(test); 
   t= norm(test);
   t(t<1)=1;      
   Dictionary(:,k+1)=test/t;
        
  
end 



dcost=zeros(1,20);
ccost=zeros(1,20);
rcost=zeros(1,20);
aPsnr=zeros(1,20);

allcost=zeros(1,20);
rho = 0.05;
tau = 0.001;
alpha=zeros(K,p);
New=I;
for i=1:20
      
         X=Expatch(n,New);
        
         alpha=fista(alpha,X,Dictionary,rho,tau);
         
         
         Dictionary=fistaD(X,alpha, Dictionary,rho);
         New = fistaX(New,I,Dictionary,alpha,rho);
    
    
    fprintf("calculating"+ i);
    dcost(1,i)=sum(vecnorm(I-New).^2);
    ccost(1,i)=rho*sum(vecnorm(X-Dictionary*alpha).^2);
    rcost(1,i)=tau*norm(alpha(:),1);
    aPsnr(1,i)=psnr(New,Io);
    allcost(1,i)=dcost(1,i)+ccost(1,i)+ rcost(1,i);
    
    %%% track PSNR
%     outs.snr(indIter) = evaluateSNR(x, imcrop(xhat,position));
    
    %%% track cost
    %%% plot global state
    figure(101);
    set(gcf, 'Color', 'w');
    
    subplot(3, 4, [1 2 5 6]);
    drawKernels(Dictionary);
    axis equal off tight;
    colormap gray;
    title('Dictionary');
    set(gca, 'FontSize', 18);
    
    subplot(3, 4, 3);
    imagesc(New);
    axis equal off tight;
    colormap gray;
    title(sprintf('xhat: %.2f dB', psnr(New, Io)));
    set(gca, 'FontSize', 18);
    
    subplot(3, 4, 4);
    imagesc(I);
    axis equal off tight;
    colormap gray;
    title(sprintf('xhat: %.2f dB', psnr(I, Io)));
    set(gca, 'FontSize', 18);
    
    subplot(3, 4, 7:8);
   
    semilogy(1:i, aPsnr(1:i), 'b-', 'LineWidth', 1.5);
    xlim([1 20]);
    grid on;
    title('psnr');
    set(gca, 'FontSize', 15);
   
    
    subplot(3, 4, 9);
%     title(sprintf('||X-Y||^2'));
    semilogy(1:i, dcost(1:i), 'b-', 'LineWidth', 1.5);
    xlim([1 20]);
    grid on;
    title('||X-Y||^2');
    set(gca, 'FontSize', 12);
    
    subplot(3, 4, 10);
%     title(sprintf("||RX-Dalpha||"));
    semilogy(1:i, ccost(1:i), 'c-', 'LineWidth', 1.5);
    xlim([1 20]);
    grid on;
     title('||RX-Dalpha||^2_2');
    set(gca, 'FontSize', 12);
    
    subplot(3, 4, 11);
%     title(sprintf('|alpha|_1'));
    semilogy(1:i, rcost(1:i), 'm-', 'LineWidth', 1.5);
    xlim([1 20]);
    grid on;
     title('|alpha|_1');
    set(gca, 'FontSize', 12);
    
     subplot(3, 4, 12);
%     title(sprintf('allcost'));
    semilogy(1:i, allcost(1:i), 'm-', 'LineWidth', 1.5);
    xlim([1 20]);
    grid on;
     title('allcost');
    set(gca, 'FontSize', 12);
    

    
    drawnow;
end




