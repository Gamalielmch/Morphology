clear
padi='C:\Users\Usuario\Documents\MATLAB\particulas';
[filenamed, pathd]= uigetfile({'*.*'},'Images of CLASS ','MultiSelect', 'on',padi);
if pathd==0
    error('No files, try again, run function again');
end

if iscell(filenamed)
    nclass=length(filenamed);
    images=cell(nclass,1);
    for i=1:nclass
        I=imread([pathd, filenamed{i}]);
        if size(I,3)>3
            I=I(:,:,1:3);
        end
        images{i}=im2bw(mat2gray(I));
        
    end
else
    I=imread([pathd, filenamed]);
    if size(I,3)>3
        I=I(:,:,1:3);
    end
    images{1}=im2bw(mat2gray(I));
    nclass=1;
end
optimo=0.999976100000000;
nc=128;   %%% número de armónicos
ncr=3;    %%% Número de armónicos para la reconstrucción
graf=1;   %%% visualización de la reconstrucción
ii = 1: nc;
% figure, hold on
for i=1:nclass
    
    [a,b,c,d,A0,C0,mn]=elifou(images{i},nc,1);
    real=zeros(mn,2);
    recon=real;
    for t = 1 : mn
        real(t,1) =sum( (a(ii) .* cos(2 .* ii .* pi .* t / mn) + b(ii) .* sin(2 .* ii .* pi .* t / mn)));
        real(t,2) = sum((c(ii) .* cos(2 .* ii .* pi .* t / mn) + d(ii) .* sin(2 .* ii .* pi .* t / mn)));
    end
    
    fc=1e-10:0.01:0.9;
   % figure
    for um=1:length(fc)
        coef= fir1(200,fc(um),'low');
        [h,~] = freqz(coef,1,length(a));
        h=abs(h)';
        for t = 1 : mn
            recon(t,1) = sum(h.*(a.* cos(2 *  ii.* pi * t / mn) + b.* sin(2 * ii.* pi * t / mn) )  );
            recon(t,2) = sum(h.*(c.* cos(2 * ii.* pi * t / mn) + d.* sin(2 * ii.* pi * t / mn) )  );
        end
        [NSx, ~] = nashsutcliffe([(1:t)', recon(:,1)], [(1:t)', real(:,1)]);
        [NSy, ~] = nashsutcliffe([(1:t)', recon(:,2)], [(1:t)', real(:,2)]);
         umbral(um,i)=sqrt((NSx^2+NSy^2));
%          umbral(um, i)=sum(   sqrt(  (recon(:,1)- real(:,1)).^2 +  (recon(:,2)- real(:,2)).^2  )  ) /...
%              sum(   sqrt(  (recon(:,1)).^2 +  (recon(:,2)).^2  )  );
    end
    
end


l=-10;
optimo=0;
figure
dif=zeros(1,size(umbral,1));
for opt=1:size(umbral,1)
    mini=min(umbral(opt,1:20));
    maxi=max(umbral(opt,21:25));
    dif(opt)=mini-maxi;
    if (mini-maxi)>l
        l=(mini-maxi);
        optimo=opt;
    end
end
for i=1:20
group{i,1}='class';
end
for i=21:25
group{i,1}='well roundness';
end
svmStruct = svmtrain(data,group,'ShowPlot',true);
%class well rounded fc=0.024
%class rounded fc=0.04691
%class sub-rounded fc=0.1404
%class angular fc=0.1405
% 
% figure
% plot([1:5],umbral(optimo,1:5),'r+')
% hold on
% plot([6:10],umbral(optimo,6:10),'g+')
% plot([11:15],umbral(optimo,11:15),'b+')
% plot([16:20],umbral(optimo,16:20),'k+')
% plot([21:25],umbral(optimo,21:25),'c+')
