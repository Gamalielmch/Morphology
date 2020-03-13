%%% Frecuencias bunas: 0.0166 0.0431


%class well rounded fc=0.024
%class rounded fc=0.04691
%class sub-rounded fc=0.1404
%class angular fc=0.1404
% % 
% % clear
% % padi='C:\Users\Usuario\Documents\MATLAB\particulas';
% % [filenamed, pathd]= uigetfile({'*.*'},'Images of CLASS ','MultiSelect', 'on',padi);
% % if pathd==0
% %     error('No files, try again, run function again');
% % end
% % 
% % assi=[ones(1,5), ones(1,5)*2, ones(1,5)*3, ones(1,5)*4, ones(1,5)*5]';
% % 
% % if iscell(filenamed)
% %     nclass=length(filenamed);
% %     images=cell(nclass,1);
% %     for i=1:nclass
% %         I=imread([pathd, filenamed{i}]);
% %         if size(I,3)>3
% %             I=I(:,:,1:3);
% %         end
% %         images{i}=im2bw(mat2gray(I));
% %         
% %     end
% % else
% %     I=imread([pathd, filenamed]);
% %     if size(I,3)>3
% %         I=I(:,:,1:3);
% %     end
% %     images{1}=im2bw(mat2gray(I));
% %     nclass=1;
% % end
% % nc=128;   %%% número de armónicos
% % ncr=3;    %%% Número de armónicos para la reconstrucción
% % graf=1;   %%% visualización de la reconstrucción
% % ii = 1: nc;
% % % figure, hold on
% % a=zeros(nclass,nc);
% % b=zeros(nclass,nc);
% % c=zeros(nclass,nc);
% % d=zeros(nclass,nc);
% % mn=zeros(nc,1);
% % for i=1:nclass
% %     
% %     [a(i,:),b(i,:),c(i,:),d(i,:),~,~,mn(i)]=elifou(images{i},nc,1);
% %     realt=zeros(mn(i),2);
% %     for t = 1 : mn(i)
% %         realt(t,1) =sum( (a(i,ii)  .* cos(2 .* ii .* pi .* t / mn(i)) + b(i,ii)  .* sin(2 .* ii .* pi .* t ./ mn(i))));
% %         realt(t,2) = sum((c(i,ii)  .* cos(2 .* ii .* pi .* t / mn(i)) + d(i,ii)  .* sin(2 .* ii .* pi .* t ./ mn(i))));
% %     end
% %     %     plot(realt(:,1), realt(:,2), 'b', 'linewidth', 3);
% %     real{i}={realt};
% % end
% % 
% % 
% % 
% % fc=1e-6:0.001:0.9;
% % C= combnk(fc,2);
% % fc1=0.01:0.0001:0.040;
% % fc2=0.0290:0.0001:0.06;
% % 
% % umbral1=zeros(1, 25);
% % umbral2=zeros(1, 25);
% % med=zeros(2,5);
% % fco=zeros(2,5);
% % factible=0;
% % umb=3;
% % umb2=10;
% % % figure
% % for um=1:length(fc1)
% %     coef= fir1(200,fc1(um),'low');
% %     [h,~] = freqz(coef,1,length(a));
% %     h1=abs(h)';
% %     for um2=1:length(fc2)
% %         coef= fir1(200,fc2(um2),'low');
% %         [h,~] = freqz(coef,1,length(a));
% %         h2=abs(h)';
% %         
% %         for i=1:nclass
% %             recon=zeros(mn(i),2);
% %             recon2=zeros(mn(i),2);
% %             recon(:,1)=sum(h1.*(a(i,ii).* cos(2 *  ii.* pi .* [1 : mn(i)]' / mn(i)) + b(i,ii).* sin(2 * ii.* pi .* [1 : mn(i)]' / mn(i)) ),2  );
% %             recon(:,2)=sum(h1.*(c(i,ii).* cos(2 *  ii.* pi .* [1 : mn(i)]' / mn(i)) + d(i,ii).* sin(2 * ii.* pi .* [1 : mn(i)]' / mn(i)) ),2  );
% %             recon2(:,1)=sum(h2.*(a(i,ii).* cos(2 *  ii.* pi .* [1 : mn(i)]' / mn(i)) + b(i,ii).* sin(2 * ii.* pi .* [1 : mn(i)]' / mn(i)) ),2  );
% %             recon2(:,2)=sum(h2.*(c(i,ii).* cos(2 *  ii.* pi .* [1 : mn(i)]' / mn(i)) + d(i,ii).* sin(2 * ii.* pi .* [1 : mn(i)]' / mn(i)) ),2  );
% %             realt=real{i};
% %             realt=realt{1};
% %             umbral1(i)=sum(   sqrt(  (recon(:,1)- realt(:,1)).^2 +  (recon(:,2)- realt(:,2)).^2  )  ) /...
% %                 sum(   sqrt(  (recon(:,1)).^2 +  (recon(:,2)).^2  )  );
% %             umbral2(i)=sum(   sqrt(  (recon2(:,1)- realt(:,1)).^2 +  (recon2(:,2)- realt(:,2)).^2  )  ) /...
% %                 sum(   sqrt(  (recon2(:,1)).^2 +  (recon2(:,2)).^2  )  );
% %         end
% %         
% %         
% %         
% %         oo=1;
% %         for ij=0:5:20
% %             med(1,oo)=mean(umbral1(ij+1:ij+5));
% %             med(2,oo)=mean(umbral2(ij+1:ij+5));
% %             oo=oo+1;
% %         end
% %         idx = kmeans([umbral1;umbral2]',5, 'Start', med');
% %         if sum(abs(idx-assi))<umb
% %             factible=factible+1;
% %             umb2(factible)=sum(abs(idx-assi));
% %             fco(:,factible)=[fc1(um),fc2(um2)];
% %         end
% %     end
% %     
% % end
% % save ('result3.mat')



opt=-200;
longi=length(fco);
fprintf('%5d iteration de',longi) 
for j=1:length(fco)
coef= fir1(200,fco(1,j),'low');
[h,~] = freqz(coef,1,length(a));
h1=abs(h)';

coef= fir1(200,fco(2,j),'low');
[h,~] = freqz(coef,1,length(a));
h2=abs(h)';

for i=1:nclass
    recon=zeros(mn(i),2);
    recon2=zeros(mn(i),2);
    recon(:,1)=sum(h1.*(a(i,ii).* cos(2 *  ii.* pi .* [1 : mn(i)]' / mn(i)) + b(i,ii).* sin(2 * ii.* pi .* [1 : mn(i)]' / mn(i)) ),2  );
    recon(:,2)=sum(h1.*(c(i,ii).* cos(2 *  ii.* pi .* [1 : mn(i)]' / mn(i)) + d(i,ii).* sin(2 * ii.* pi .* [1 : mn(i)]' / mn(i)) ),2  );
    recon2(:,1)=sum(h2.*(a(i,ii).* cos(2 *  ii.* pi .* [1 : mn(i)]' / mn(i)) + b(i,ii).* sin(2 * ii.* pi .* [1 : mn(i)]' / mn(i)) ),2  );
    recon2(:,2)=sum(h2.*(c(i,ii).* cos(2 *  ii.* pi .* [1 : mn(i)]' / mn(i)) + d(i,ii).* sin(2 * ii.* pi .* [1 : mn(i)]' / mn(i)) ),2  );
    realt=real{i};
    realt=realt{1};
    umbral1(i)=sum(   sqrt(  (recon(:,1)- realt(:,1)).^2 +  (recon(:,2)- realt(:,2)).^2  )  ) /...
        sum(   sqrt(  (recon(:,1)).^2 +  (recon(:,2)).^2  )  );
    umbral2(i)=sum(   sqrt(  (recon2(:,1)- realt(:,1)).^2 +  (recon2(:,2)- realt(:,2)).^2  )  ) /...
        sum(   sqrt(  (recon2(:,1)).^2 +  (recon2(:,2)).^2  )  );
end

mini=min(umbral2(1:10));
maxi=max(umbral2(11:15));
if (mini-maxi)> opt
    opt=(mini-maxi);
    bueno=j;
end
fprintf('\b\b\b\b%4.0f%',j)
drawnow
% 
% figure
% plot(umbral1(1:5),umbral2(1:5),'r+')
% hold on
% plot(umbral1(6:10),umbral2(6:10),'g+')
% plot(umbral1(11:15),umbral2(11:15),'b+')
% plot(umbral1(16:20),umbral2(16:20),'k+')
% plot(umbral1(21:25),umbral2(21:25),'c+')
% close
end







% hold on
% th = 0:pi/50:2*pi;
% xunit = r * cos(th) + x;
% yunit = r * sin(th) + y;
% h = plot(xunit, yunit);
% hold off
%
%
%
%
%
%
%
%
%
% l=-10;
% optimo=0;
% figure
% dif=zeros(1,size(umbral,1));
% for opt=1:size(umbral,1)
%     mini=min(umbral(opt,1:20));
%     maxi=max(umbral(opt,21:25));
%     dif(opt)=mini-maxi;
%     if (mini-maxi)>l
%         l=(mini-maxi);
%         optimo=opt;
%     end
% end
%
%
% figure
% plot([1:5],umbral(optimo,1:5),'r+')
% hold on
% plot([6:10],umbral(optimo,6:10),'g+')
% plot([11:15],umbral(optimo,11:15),'b+')
% plot([16:20],umbral(optimo,16:20),'k+')
% plot([21:25],umbral(optimo,21:25),'c+')
%
%
%
% [pks,locs] = findpeaks(dif,'SortStr','descend');
% %class well rounded fc=0.024
% %class rounded fc=0.04691
% %class sub-rounded fc=0.1404
% %class angular fc=0.1405
% plot(dif), hold on
% plot(locs,pks,'ro')
% text(locs+.02,pks,num2str((1:numel(pks))'))
%
% data(1,:)=umbral(locs(1),:);
% data(2,:)=umbral(locs(2),:);
%
% for i=1:20
%     group{i,1}='class';
% end
% for i=21:25
%     group{i,1}='well roundness';
% end
% figure,
% svmStruct = svmtrain(data,group,'ShowPlot',true);








