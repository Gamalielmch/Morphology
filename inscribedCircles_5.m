function roundness = inscribedCircles_5 (img_route,armon,dist_max,im_name,folder,display)
%display=0 no se muestra nada, display=1 se muestra contorno, diaplay=2 se
%muestra contornos e imagen
%close all;
% read and binarization image

%%%%%%%%% Read and area normalization 
im = imread(img_route);
im = imbinarize(im(:,:,1));
im=regionprops(im,'Image');
im=im(1).Image;
fac=450/size(im,1);
im=imresize(im,fac);
im = padarray(im,[15 15],'both');


%%%%%%%%% spur contour remove 
windowSize = 5;
kernel = ones(windowSize) / windowSize ^ 2;
im = conv2(single(im), kernel, 'same');
im = im > 0.5; % Rethreshold

%%%% rotate 
BW = im;
BW=imrotate(BW,90);
[nx,ny]=size(BW);

%Finding maximum circunscribed circle
[cx,cy,radii]=max_circun_circle(BW,0);

% Fourier eliptico smoothing
imSmooth = bw_chain_fil(BW,armon,0);

x=imSmooth(:,2); y=imSmooth(:,1);
temp=false(nx,ny);
L=round(x);
L(:,2)=round(y);
for i=1:length(x)
    temp(L(i,1),L(i,2))=1;
end

temp=imfill(bwmorph(temp,'bridge'),'holes');

X = [x,y];
[~,~,K2] = curvature(X);
% mag=sqrt(K2(:,1).^2+K2(:,2).^2);
Si=sign(K2);
temp2=false(size(temp));
LX=false(length(x),1);
xt=[]; yt=[];
points=[];
for i=2:length(x)-1
    temp2(L(i,1)+Si(i,1),L(i,2)+Si(i,2))=1;
    t=and(temp2,temp);
    LX(i)=sum(t(:));
    if LX(i)>0 %&& mag(i)>0.0008
        xt=[xt x(i)];
        yt=[yt y(i)];
        points=[points, i];
    end
    temp2(L(i,1)+Si(i,1),L(i,2)+Si(i,2))=0;
end



[results]=fitting_circ(xt,yt,radii,dist_max,imSmooth,0);
if size(results,1)>7
[~,ind]= min([results{:,1}]);
results(ind,:)=[];
[~,ind]= max([results{:,1}]);
results(ind,:)=[];
[~,ind]= min([results{:,1}]);
results(ind,:)=[];
[~,ind]= max([results{:,1}]);
results(ind,:)=[];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% show results
if display==0
    
        roundness = mean([results{:,1}])/radii;
        if roundness>1
            roundness=1;
        end
        return
end
if display>1
fig1=figure('color',[1 1 1]);
axes1 = axes('Parent',fig1);
% BW = im;
imshow(flip(~im,2))
hold on
end
%%%%% contour
fig2=figure('color',[1 1 1]);
axes2 = axes('Parent',fig2);
plot(imSmooth(:,2),imSmooth(:,1),'k')
hold on
% selected points
plot(imSmooth(points,2),imSmooth(points,1),'r.');
% plot maximum circunscribed circle
theta=linspace(0,2*pi,2000);
[xr,yr]=pol2cart(theta, repmat(radii,1,2000));
plot(yr+cx,xr+cy,'b')
if display>1
plot(axes1,yr+cx,xr+cy,'y')
end
% plot curvature circle
for i=1:size(results,1)
    R=results{i,1};
    if(~isempty(R))
        cent=results{i,2};
        rt=R*ones(1,2000);
        [xtt,yyt]=pol2cart(theta,rt);
        figure(fig2)
        plot(axes2,cent(1),cent(2),'r.')
        plot(axes2,xtt+cent(1),yyt+cent(2),'color',[0.5,0.25,0])
        if display>1
        plot(axes1,xtt+cent(1),yyt+cent(2),'g')
        end
    end
end


% dim = [.15 .85 .05 .05];
dim=[0.125 0.923809525015809 0.253571421546596 0.0642857130794299];
roundness = mean([results{:,1}])/radii;
if(roundness > 1)
   roundness = 1; 
end
view(axes2,[180 90]);
if display>1
view(axes1,[-180 -90]);
end
annotation('textbox',dim,'String',"Roundness:"+roundness,'FitBoxToText','on');
axis equal
saveas(fig2,[folder,'/',im_name])
% close (fig1)
% close (fig2)
end


%maximum circunscribed circle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [cx,cy,radii]=max_circun_circle(BW,show)
dist = bwdist(~BW);
radii=max(dist(:));
[cx,cy]=find(dist==radii);
cx=cx(1);
cy=cy(1);
if show==1
    figure, imshow(BW), hold on
    [xr,yr]=pol2cart(linspace(0,2*pi,2000), repmat(radii,1,2000));
    plot(cy,cx,'r*', 'MarkerSize',20)
    plot(xr+cy,yr+cx)
end

end


%Distance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pd=distance(x,y,maxd,pu)
pd=[];
for i=1:length(pu)
    pin=[x(pu(i)),y(pu(i))];
    d=sqrt((pin(1)-x).^2 + (pin(2)-y).^2);
    d=find(d<maxd);
    pd=[pd,d];
end
end

%fitting_circule
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [result]=fitting_circ(x,y,radii,dist_max,imSmooth,show)

if show
    figure('color',[1 1 1]);
    plot(imSmooth(:,2),imSmooth(:,1),'k')
    hold on
    plot(x,y,'rs');
    theta=linspace(0,2*pi,2000);
end
result=cell(1,2);
cont=1;
lenx=length(x);
while lenx>2
    
    re=1;
    dif=1;
    d2=1;
    while dif>0
        re=distance(x,y,dist_max,re);
        re=unique(re);
        dif=length(re)-length(d2);
        d2=re;
    end
    
    if length(re)>15
        
        [xc,yc,R,~] = circfit(x(re),y(re));
        if R<=2*radii
            
            [~,gof,~]=fit(x(re)',y(re)','poly1');
            
            if gof.rsquare<0.94
                cent=[xc,yc];
                result(cont,1)={R};
                result(cont,2)={cent};
                cont=cont+1;
                
                if show
                    plot(x(re),y(re),'bs');
                    plot(cent(1),cent(2),'r*')
                    rt=R*ones(1,2000);
                    [xtt,yyt]=pol2cart(theta,rt);
                    plot(xtt+cent(1),yyt+cent(2),'color',[0.5,0.25,0])
                end
            end
        end
    end
    x(re)=[];
    y(re)=[];
    lenx=length(x);
end


end

%finding the best circule
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function   [xc,yc,R,a] = circfit(x,y)
x=x(:); y=y(:);
a=[x y ones(size(x))]\[-(x.^2+y.^2)];
xc = -.5*a(1);
yc = -.5*a(2);
R  =  sqrt((a(1)^2+a(2)^2)/4-a(3));
end


% function   indmax = puntos_lejanos(x,y)
% dist=zeros(length(x));
% maximo=0;
% indmax=[];
% for i=1:length(x)
%     dist(:,i)= sqrt((x(i)-x).^2 + (y(i)-y).^2);
%     if max(dist(:,i))>maximo
%         [maximo,im]=max(dist(:,i));
%         indmax=[i,im];
%     end
% end
% end

