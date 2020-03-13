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

nc=128;   %%% número de armónicos
ncr=3;    %%% Número de armónicos para la reconstrucción
graf=1;   %%% visualización de la reconstrucción
% figure, hold on
for i=1:nclass
        image=images{i};
        imt=uint8(image*255);
        s = regionprops(image,'centroid','PixelIdxList');
        cen = cat(1,s.Centroid);
        obj=length(s);
    for io=1:obj
        image=false(size(image));
        image(s(io).PixelIdxList)=1;
        [a,b,c,d,A0,C0,mn]=elifou(image,nc,1);
        real=zeros(mn,2);
        recon=real;
        mag=zeros(1, mn);
        for nt=1:nc
            thn=0.5*atand ( (2*(a(nt)*b(nt) +c(nt)*d(nt))) / (a(nt).^2+c(nt).^2-b(nt).^2-d(nt).^2)  );
            t=[cosd(thn) sind(thn); -sind(thn) cosd(thn)]*[a(nt),c(nt);b(nt),d(nt)];
            ang1=atand(t(3)/t(1));
            ang2=atand(t(4)/t(2));
            l1=sqrt(t(1).^2 + t(3).^2);
            l2=sqrt(t(2).^2 + t(4).^2);
            area(nt)= abs(l1) * abs(l2)  *pi;
            [ma,ma2]=max(abs([l1, l2]));
            if ma2==1
                exc(nt)= (l2/l1);
            else
                exc(nt)= (l1/l2);
            end
        end
        for ip=4:length(area)-5
            suma(ip-3)=sum(area(ip:ip+5));
        end
        area=area/sum(area);
        %    figure,% subplot(1,2,1)
        areacum=(cumsum(area));
        %    plot(log(1:length(suma)),log(suma))
        um=find(areacum>=0.993,1);
        
        %    p=plot(log(1:length(areacum)),log(areacum));
        %    [ind,ix]=knee_pt(log(1:length(areacum)),log(areacum));
        %    sum(exc(ind+1:end).*area(ind+1:end))
        %       figure
        %       plot(areacum(1:15))
        %      hold on
        %      plot(um,areacum(um),'o');
        %      title(num2str(um))
        %     %
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
            for ii = 1: um
                x = x + (a(ii) * cos(2 * ii * pi * t / mn) + b(ii) * sin(2 * ii * pi * t / mn));
                y = y + (c(ii) * cos(2 * ii * pi * t / mn) + d(ii) * sin(2 * ii * pi * t / mn));
            end
            recon(t,1) = A0 + x;
            recon(t,2) = C0 + y;
            mag(t)=sqrt((real(t,1)-recon(t,1))^2+(real(t,2)-recon(t,2))^2);
            mag(t)=mag(t)*sign(real(t,2)-recon(t,2));
        end

        figure
        plot(recon(:,1), recon(:,2), 'r', 'linewidth', 2);
        axis equal
        hold on
        plot(real(:,1), real(:,2), 'k', 'linewidth', 2);
        a = abs(polyarea(recon(:,1),recon(:,2))-polyarea(real(:,1),real(:,2)))/(polyarea(real(:,1),real(:,2)))*100;
        imt=insertText(imt,[cen(io,1),cen(io,2)],num2str(a),'FontSize',45,'BoxColor',...
                       'yellow','BoxOpacity',0.4,'TextColor','black');
%         title(num2str(a))
    end
     figure, imshow(imt);
end

%   figure, subplot(1,2,1)
%      axis equal
%      hold on
%      plot(real(:,1), real(:,2), 'k', 'linewidth', 2);
%     for t = 1 : mn
%          plot(recon(t,1), recon(t,2), 'r.', 'linewidth', 3);
%     end
