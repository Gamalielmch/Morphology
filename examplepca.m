%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%           REDUCE EL NUMERO DE SERIES                %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%                 Carga de imágenes                       %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        images{i}=im2bw(I);
        
    end
else
    I=imread([pathd, filenamed]);
    %     if size(I,3)>3
    %         I=I(:,:,1:3);
    %     end
    images{1}=im2bw(I);
    nclass=1;
end

npc=2; %%% número de componentes principales
nc=339; %%% número de armónicos
corte=zeros(1,nclass);
graf=1;
figure
for i=1:nclass
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%     Calculo de los coeficientes de Fourier eliptico     %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [a,b,c,d,A0,C0,mn]=elifou(images{i},nc,1);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Reconstrucción con todos los coeficientes sin análisis PCA %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if graf
        recon=zeros(mn,2);
        for t = 1 : mn
            x = 0.0;
            y = 0.0;
            for ii = 1: nc
                x = x + (a(ii) * cos(2 * ii * pi * t / mn) + b(ii) * sin(2 * ii * pi * t / mn));
                y = y + (c(ii) * cos(2 * ii * pi * t / mn) + d(ii) * sin(2 * ii * pi * t / mn));
            end
            recon(t,1) = A0 + x;
            recon(t,2) = C0 + y;
        end
        obsData=[(1 : mn)',recon];
        figr = figure;
        plot(recon(:,1), recon(:,2), 'b', 'linewidth', 3);
        axis equal
        hold on
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%                    Análisis PCA                            %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    dat=[a;b;c;d];
    [r,c1]=size(dat);
    m=mean(dat,2);
    data=dat-repmat(m,1,c1);
    cova=cov(data');
    [eigvector,eigvl]=eig(cova);
    eigvl=diag(eigvl);
    eigvl=eigvl(end:-1:1);
    %     sum(eigvl(1:npc))/(sum(eigvl))
    
    Data_xD=eigvector(:,end-npc+1:end)'*data;
    [vy,vx]=ecdf(abs(Data_xD(2:30)));
    plot(vx,vy)
    hold on
    % % npc dimensional case
    Res_xD= (eigvector(:,end-npc+1:end))*Data_xD;
    Res_xD=Res_xD+repmat(m,1,c1);
    %plot(Data_xD(2:50))
    a=Res_xD(1,:);
    b=Res_xD(2,:);
    c=Res_xD(3,:);
    d=Res_xD(4,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Reconstrucción con todos los coeficientes con análisis PCA %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if graf
        
        recon=zeros(mn,2);
        for t = 1 : mn
            x = 0.0;
            y = 0.0;
            for ii = 1: nc
                x = x + (a(ii) * cos(2 * ii * pi * t / mn) + b(ii) * sin(2 * ii * pi * t / mn));
                y = y + (c(ii) * cos(2 * ii * pi * t / mn) + d(ii) * sin(2 * ii * pi * t / mn));
            end
            recon(t,1) = A0 + x;
            recon(t,2) = C0 + y;
        end
        simData=[(1 : mn)',recon];
        recon = [recon; recon(1,1) recon(1,2)];
        figure(figr), plot(recon(:,1), recon(:,2), 'r', 'linewidth', 2);
        [NSout(i), ~] = nashsutcliffe(obsData, simData);
        
    end
end
