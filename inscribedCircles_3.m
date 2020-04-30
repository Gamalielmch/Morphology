function inscribedCircles_3
close all;
% read and binarization image
im = imread('testimage/d1.jpg');
im=imresize(im,1);
BW = imbinarize(im(:,:,1));

%Finding maximum circunscribed circle
[cx,cy,radii]=max_circun_circle(BW,0);

% Fourier eliptico smoothing
imSmooth = bw_chain_fil(BW,60,0);

% Getting points of curvature
points = points_curvature(imSmooth,0.20,0);

% removing convex curvature points
Seleted_points=points_convex(imSmooth,points,cx,cy,0);
points = points(Seleted_points); %Se eliminan los points que no cumplen
x=imSmooth(points,2); y=imSmooth(points,1);

% fitting circles
[results]=fitting_circ(x,y,radii,25);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% show results
figure
%%%%% contour
plot(imSmooth(:,2),imSmooth(:,1),'k')
hold on
% selected points
plot(imSmooth(points,2),imSmooth(points,1),'rx');
% plot maximum circunscribed circle
theta=linspace(0,2*pi,2000);
[xr,yr]=pol2cart(theta, repmat(radii,1,2000));
plot(yr+cx,xr+cy,'b')
% plot curvature circle
for i=1:size(results,1)
    R=results{i,1};
    cent=results{i,2};
    rt=R*ones(1,2000);
    [xtt,yyt]=pol2cart(theta,rt);
    plot(cent(1),cent(2),'r*')
    plot([cent(1),cent(1)+R],[cent(2),cent(2)],'b')
    d=sqrt((imSmooth(:,2)-cent(1)).^2 + (imSmooth(:,1)-cent(2)).^2 );
    [d,minr]=min(d);
    plot([cent(1),imSmooth(minr,2)],[cent(2),imSmooth(minr,1)],'r')
    plot(xtt+cent(1),yyt+cent(2),'g--')
end


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

%Getting curvature points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function points=points_curvature(imSmooth,delta,show)
len = size(imSmooth,1); %numero de pixeles en el contorno
points = 1;
i = 1;
cont = 1;
while i<len
    for xi = i+1:len
        c = polyfit([i,xi],[imSmooth(i,2),imSmooth(xi,2)],1);
        ylin = polyval(c,i:xi)';
        maxi_div=max(abs(imSmooth(i:xi,2)-ylin));
        if maxi_div>delta
            cont=cont+1;
            points(cont)=xi;
            i = xi;
            break
        end
    end
    i=i+1;
end
if show==1
    figure
    plot(imSmooth(:,2),imSmooth(:,1),'k')
    hold on
    for i=1:length(points)
        plot(imSmooth(points(i),2),imSmooth(points(i),1),'bs')
    end
    for  i = 1:length(points)-1
        plot([imSmooth(points(i),2),imSmooth(points(i+1),2)],[imSmooth(points(i),1),imSmooth(points(i+1),1)],'r')
    end
    plot([imSmooth(points(1),2),imSmooth(points(end),2)],[imSmooth(points(1),1),imSmooth(points(end),1)],'r')
    
end
end


% Seleted_points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Seleted_points=points_convex(imSmooth,points,cx,cy,show)

cont1=1;
for  i = 1:length(points)
    if i==1
        c=imSmooth(points(end),:);
        e=imSmooth(points(i),:);
        d=imSmooth(points(i+1),:);
    elseif i==length(points)
        c=imSmooth(points(i-1),:);
        e=imSmooth(points(i),:);
        d=imSmooth(points(1),:);
    else
        c=imSmooth(points(i-1),:);
        e=imSmooth(points(i),:);
        d=imSmooth(points(i+1),:);
    end
    
    coefficients = polyfit([cx, e(2)], [cy, e(1)], 1);
    xt=linspace(cx, e(2),1000);
    a = coefficients (1);
    b = coefficients (2);
    yt=a*xt+b;
    coefficients = polyfit([c(2), d(2)], [c(1), d(1)], 1);
    xcd=linspace(c(2), d(2),1000);
    a = coefficients (1);
    b = coefficients (2);
    ycd=a*xcd+b;
    
    [xi,~] = polyxpoly(xt,yt,xcd,ycd);
    if ~isempty(xi)
        Seleted_points(cont1)=i;
        cont1=cont1+1;
        color=[0 1 0];
    else
        
        color=[1 0 0];
    end
    if show
        figure
        plot(imSmooth(:,2),imSmooth(:,1),'k')
        hold on
        for ii=1:length(points)
            plot(imSmooth(points(ii),2),imSmooth(points(ii),1),'bs')
        end
        plot([c(2),e(2)],[c(1),e(1)],'color',color)
        plot([e(2),d(2)],[e(1),d(1)],'color',color)
        plot([c(2),d(2)],[c(1),d(1)],'color',color)
        plot(cx,cy,'r*', 'MarkerSize',20)
        plot([cx,e(2)],[cy,e(1)],'k')
        pause(2)
        close
    end
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
    pd=[pd,d'];
end
end

%fitting_circule
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [result]=fitting_circ(x,y,radii,dist_max)

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
    if length(re)>2
        [xc,yc,R,~] = circfit(x(re),y(re));
        if R<=radii
            cent=[xc,yc];
            result(cont,1)={R};
            result(cont,2)={cent};
            cont=cont+1;
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