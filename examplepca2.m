%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%           REDUCE EL NUMERO DE ARMÓNICOS                 %%%%%%
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

npc=128;    %%% número de componentes principales
nc=128;    %%% número de armónicos
graf=1;   %%% visualización de la reconstrucción
Data_xD=zeros(npc,4,nclass);
% figure, hold on
enea=zeros(nclass,80,2);
for i=1:nclass
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%     Calculo de los coeficientes de Fourier eliptico     %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [a,b,c,d,A0,C0,mn]=elifou(images{i},nc,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Reconstrucción con todos los coeficientes sin análisis PCA %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     nt=1;
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
    recon = [recon; recon(1,1) recon(1,2)];
    plot(recon(:,1), recon(:,2), 'b', 'linewidth', 3);
    axis equal
    hold on
    
   recon=zeros(mn,2);
   for t = 1 : mn
            x = 0.0;
            y = 0.0;
            for ii = 1: 40
                x = x + (a(ii) * cos(2 * ii * pi * t / mn) + b(ii) * sin(2 * ii * pi * t / mn));
                y = y + (c(ii) * cos(2 * ii * pi * t / mn) + d(ii) * sin(2 * ii * pi * t / mn));
            end
            recon(t,1) = A0 + x;
            recon(t,2) = C0 + y;
    end
    recon = [recon; recon(1,1) recon(1,2)];
    plot(recon(:,1), recon(:,2), 'r', 'linewidth', 3);
    
    if graf
        for nt=1:50
            %     recon=zeros(mn,2);
            %    for t = 1 : mn
            %                 x = 0.0;
            %                 y = 0.0;
            %                 for ii = nt
            %                     x = x + (a(ii) * cos(2 * ii * pi * t / mn) + b(ii) * sin(2 * ii * pi * t / mn));
            %                     y = y + (c(ii) * cos(2 * ii * pi * t / mn) + d(ii) * sin(2 * ii * pi * t / mn));
            %                 end
            %                 recon(t,1) = A0 + x;
            %                 recon(t,2) = C0 + y;
            %             end
            %             recon = [recon; recon(1,1) recon(1,2)];
            %  obsData=[(1 : mn)',recon];
            %         recon = [recon; recon(1,1) recon(1,2)];
            %         figr = figure;
            %              plot(recon(:,1), recon(:,2), 'b', 'linewidth', 3);
            %              axis equal
            %             hold on
            
            thn=0.5*atand ( (2*(a(nt)*b(nt) +c(nt)*d(nt))) / (a(nt).^2+c(nt).^2-b(nt).^2-d(nt).^2)  );
            t=[cosd(thn) sind(thn); -sind(thn) cosd(thn)]*[a(nt),c(nt);b(nt),d(nt)];
            ang1=atand(t(3)/t(1));
            ang2=atand(t(4)/t(2));
            l1=sqrt(t(1).^2 + t(3).^2);
            l2=sqrt(t(2).^2 + t(4).^2);
            mag= abs(l1) * abs(l1)  *pi;
            [ma,ma2]=max(abs([l1, l2]));
            
            %             plot([0 l1*cosd(ang1)],[0 l1*sind(ang1)])
            %             plot([0 l2*cosd(ang2)],[0 l2*sind(ang2)])
            %             title(num2str(ang1))
            %             hold off
            
            if ma2==1
                exc= (l2/l1);
            else
                exc= (l1/l2);
            end
            enea(i,nt,1)=exc;
            enea(i,nt,2)=mag;
        end
    end
    bin=linspace(0,1,10);
    hb=zeros(1,length(bin));
    for jj=6:size(enea,2)
        dist    = abs(enea(i,jj,1) - bin);
        minDist = min(dist);
        idx     = find(dist == minDist);
        hb(idx) = hb(idx)+enea(i,jj,2);
    end
   % figure, bar(bin,hb)
    s(i)=sum(hb(1:3));
    plot((sort (enea(i,:,1))),'k'), hold on
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     %%%%%%                    Análisis PCA                            %%%%%
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %
    %     dat=[a;b;c;d];
    %     dat=dat';
    %     [r,c1]=size(dat);
    %     m=mean(dat,2);
    %     data=dat-repmat(m,1,c1);
    %     cova=cov(data');
    %     [eigvector,eigvl]=eig(cova);
    %     [W,pc,lar]=pca(dat','NumComponents',10);
    %     eigvl=diag(eigvl);
    %     eigvl=eigvl(end:-1:1);
    %     %eigvl(1)/(sum(eigvl))
    %     Data_xD(:,:,i)=eigvector(:,end-npc+1:end)'*data;
    %     %%%%%
    %     %eigvector(:,end-npc+1:end)'
    %     Res_xD= (eigvector(:,end-npc+1:end))*Data_xD(:,:,i);
    %     Res_xD=Res_xD+repmat(m,1,c1);
    %     a=Res_xD(:,1);
    %     b=Res_xD(:,2);
    %     c=Res_xD(:,3);
    %     d=Res_xD(:,4);
    %     %plot((squeeze(Data_xD(1,:,i)))), hold on
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     %%%%%% Reconstrucción con todos los coeficientes con análisis PCA %%%%%
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     if graf
    %         recon=zeros(mn,2);
    %         for po=1:nc
    %             for t = 1 : mn
    %                 x = 0.0;
    %                 y = 0.0;
    %                 for ii = 1: po
    %                     x = x + (a(ii) * cos(2 * ii * pi * t / mn) + b(ii) * sin(2 * ii * pi * t / mn));
    %                     y = y + (c(ii) * cos(2 * ii * pi * t / mn) + d(ii) * sin(2 * ii * pi * t / mn));
    %                 end
    %                 recon(t,1) = A0 + x;
    %                 recon(t,2) = C0 + y;
    %             end
    %             simData=[(1 : mn)',recon];
    %             %         recon = [recon; recon(1,1) recon(1,2)];
    %             %         figure(figr), plot(recon(:,1), recon(:,2), 'r', 'linewidth', 2);
    %             %         hold off
    %             [NSout(i,po), ~] = nashsutcliffe(obsData, simData);
    %         end
    %         plot(NSout(i,2:20))
    %     end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%                    Gráfica de los PCA                      %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% color=lines(5);
% figure, hold on
% Data_xD=flipdim(Data_xD,3);
% marcador=['o','+','*','x'];
% %  for j=1:5
% %      for i=1:5
% %          for n=2
% %          scatter(Data_xD(2,n,i+((j-1)*5)),Data_xD(3,n,i+((j-1)*5)),'MarkerEdgeColor',color(j,:),'Marker',marcador(n))
% %          end
% %      end
% %      j
% %  end
% %,'MarkerFaceColor',color(j,:);
% for j=1:5
%     for i=1:5
%         for n=3
%             scatter3(Data_xD(1,n,i+((j-1)*5)),Data_xD(2,n,i+((j-1)*5)),Data_xD(3,n,i+((j-1)*5)),'MarkerEdgeColor',color(j,:),'Marker',marcador(n))
%         end
%     end
% end
%
% % for j=1:1
% %     for i=1:5
% %         scatter3(Data_xD(1,:,i+((j-1)*5)),Data_xD(2,:,i+((j-1)*5)),Data_xD(3,:,i+((j-1)*5)),'MarkerFaceColor',color(i,:),'MarkerEdgeColor',color(i,:))
% %     end
% % end
% % view(gca,[-35.6 34.8]);
% % grid(gca,'on');
% %





















