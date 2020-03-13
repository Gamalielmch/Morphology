function roundness_using_curvature
close all
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

I=imresize(im2bw(images{1}),1);


perifil=bw_chain_fil(I,50,0);

figure
peri=bwboundaries(I);
peri=peri{:};
bwd=bwdist(~I);
imshow(bwd,[],'InitialMagnification','fit')
hold on 
[r,m]=max(bwd(:));
[y,x]=ind2sub(size(bwd),m);
plot(peri(:,2),peri(:,1),'r')
plot(perifil(:,1),perifil(:,2),'g')
plot(x,y,'ro','MarkerEdgeColor','b',...
    'MarkerFaceColor',[0 0 1])
dth=linspace(1e-10,2*pi,10000);
xc=x+r*cos(dth);
yc=y+r*sin(dth);
plot(xc,yc,'y')

%%% Set up axis
xmin=min(peri(:,2));
xmax=max(peri(:,2));
dx=round((xmax-xmin)*.05);
ymin=min(peri(:,1));
ymax=max(peri(:,1));
dy=round((ymax-ymin)*.05);
xlim([xmin-dx,xmax+dx])
ylim([ymin-dy,ymax+dy])
axis equal
axis ij

%%% key points
c=regionprops(I,'centroid');
c=c.Centroid;
PCD= max(sqrt((peri(:,2)-c(1)).^2+(peri(:,1)-c(2)).^2));
xc=c(1)+PCD*cos(dth);
yc=c(2)+PCD*sin(dth);
% plot(c(1),c(2),'go','MarkerEdgeColor','g',...
%     'MarkerFaceColor',[0 1 0])
plot(xc,yc,'g')
delta=PCD*0.0003;
delta
