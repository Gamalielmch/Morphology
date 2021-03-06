function roundness = inscribedCircles_4 (img_route,armon,dist_max,im_name)
close all;
% read and binarization image
    inscribedRadii = 0;
    im = imread(img_route);
    im=imresize(im,1.5);
    BW = imbinarize(im(:,:,1));
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
    imshow(temp)

    X = [x,y];
    [L2,R2,K2] = curvature(X);

    figure;
    plot(L2,R2)
    title('Curvature radius vs. cumulative curve length')
    xlabel L
    ylabel R
    figure;
    h = plot(x,y); grid on; axis equal
    set(h,'marker','.');
    xlabel x
    ylabel y
    title('2D curve with curvature vectors')
    hold on

    Si=sign(K2);
    temp2=false(size(temp));
    LX=false(length(x),1);
    xt=[]; yt=[];
    points=[];
    for i=2:length(x)-1
       temp2(L(i,1)+Si(i,1),L(i,2)+Si(i,2))=1;
       t=and(temp2,temp);
       LX(i)=sum(t(:));
       if LX(i)>0 
           xt=[xt x(i)];
           yt=[yt y(i)];
           points=[points, i];
       end
       temp2(L(i,1)+Si(i,1),L(i,2)+Si(i,2))=0;
    end

    [results]=fitting_circ(xt,yt,radii,dist_max,imSmooth,0);



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% show results
    fig1=figure('color',[1 1 1]);
    axes1 = axes('Parent',fig1);
    BW = imbinarize(im(:,:,1));
    imshow(flip(~BW,2))
    hold on
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
    plot(axes1,yr+cx,xr+cy,'y')

    % plot curvature circle
    for i=1:size(results,1)
        R=results{i,1};
        if(~isempty(R))
            inscribedRadii = inscribedRadii + R;
            cent=results{i,2};
            rt=R*ones(1,2000);
            [xtt,yyt]=pol2cart(theta,rt);
            figure(fig2)
            plot(axes2,cent(1),cent(2),'r.')
            %plot(axes2,[cent(1),cent(1)+R],[cent(2),cent(2)],'b')
            %d=sqrt((imSmooth(:,2)-cent(1)).^2 + (imSmooth(:,1)-cent(2)).^2 );
            %[d,minr]=min(d);
            %plot(axes2,[cent(1),imSmooth(minr,2)],[cent(2),imSmooth(minr,1)],'r')
            plot(axes2,xtt+cent(1),yyt+cent(2),'color',[0.5,0.25,0])
            plot(axes1,xtt+cent(1),yyt+cent(2),'g')
        end
    end



%     set (axes1, 'xdir', 'reverse' )
%     axes(axes2)
%     axis equal
%     view(axes2,[-180 90]);
    dim = [.15 .85 .05 .05];
    roundness = (inscribedRadii/size(results,1))/radii;
    annotation('textbox',dim,'String',"Roundness:"+round(roundness*10)/10,'FitBoxToText','on');
    saveas(fig2,"resultados/"+im_name)
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

        if length(re)>2

            [xc,yc,R,~] = circfit(x(re),y(re));
            if R<=0.8*radii
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