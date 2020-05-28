function  [spher]=sphericity(img_route,show)
    %close all;

    im = imread(img_route);
    im=imresize(im,1.5);
    BW = imbinarize(im(:,:,1));
    BW=imrotate(BW,90);
    
    dist = bwdist(~BW);
    radii=max(dist(:));
    [cx,cy]=find(dist==radii);
    cx=cx(1);
    cy=cy(1);
    imSmooth = bw_chain_fil(BW,40,0);
    circle = enclosingCircle(imSmooth);
    
    if show==1
        
        figure, plot(imSmooth(:,1),imSmooth(:,2)), hold on
        [xr,yr]=pol2cart(linspace(0,2*pi,2000), repmat(radii,1,2000));
        plot(cy,cx,'r*', 'MarkerSize',20)
        plot(xr+cy,yr+cx)
        
         [xr,yr]=pol2cart(linspace(0,2*pi,2000), repmat(circle(3),1,2000));
         plot(circle(1),circle(2),'b--o', 'MarkerSize',20)
         plot(xr+circle(1),yr+circle(2))
    end
    %spher = 1;
    spher = radii/circle(3);
end

function   [xc,yc,R,a] = circfit(x,y)
    x=x(:); y=y(:);
    a=[x y ones(size(x))]\[-(x.^2+y.^2)];
    xc = -.5*a(1);
    yc = -.5*a(2);
    R  =  sqrt((a(1)^2+a(2)^2)/4-a(3));
end