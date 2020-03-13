
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

nc=128;   %%% número de armónicos
ncr=3;    %%% Número de armónicos para la reconstrucción
graf=1;   %%% visualización de la reconstrucción
% figure, hold on
for i=1:nclass
        [a,b,c,d,A0,C0,mn]=elifou(images{i},nc,1);
        real=zeros(mn,2);
        recon=real;
        for t = 1 : mn
            x = 0.0;
            y = 0.0;
            for ii = 1: nc
                x = x + (a(ii) * cos(2 * ii * pi * t / mn) + b(ii) * sin(2 * ii * pi * t / mn));
                y = y + (c(ii) * cos(2 * ii * pi * t / mn) + d(ii) * sin(2 * ii * pi * t / mn));
            end
            real(t,1) = A0 + x;
            real(t,2) = C0 + y;
            
            x = 0.0;
            y = 0.0;
            for ii = 1: ncr
                x = x + (a(ii) * cos(2 * ii * pi * t / mn) + b(ii) * sin(2 * ii * pi * t / mn));
                y = y + (c(ii) * cos(2 * ii * pi * t / mn) + d(ii) * sin(2 * ii * pi * t / mn));
            end
            recon(t,1) = A0 + x;
            recon(t,2) = C0 + y;     
        end
        
        
        
                 figr = figure;
                 plot(real(:,1), real(:,2), 'k', 'linewidth', 3);
                 axis equal
                 hold on
                 plot(recon(:,1), recon(:,2), 'b', 'linewidth', 2);
end
