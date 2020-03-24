close all;
load franke
im = imread('im_pruebas/d4.jpg');
%im = imrotate(im,25,'bicubic');
%imshow(im);
BW = imbinarize(im(:,:,1));

%PASO 1
dist = bwdist(~BW);
radii=max(dist(:));

[x,y]=find(dist==radii);
[B,L] = bwboundaries(BW);
N = 40;
%figure,imshow(label2rgb(L, @jet, [.5 .5 .5]))
%imshow(a_bw);
for k = 1:size(B,1)
   boundary = B{k};
   %figure,plot(boundary(:,2),boundary(:,1))
end

Ibw = imfill(BW,'holes');
Ilabel = bwlabel(Ibw);
stat = regionprops(Ilabel,'centroid');
centroid = stat.Centroid;

res = boundary - [x,y]; %Se resta las posiciones del centroide, para poner la figura en el origen
[theta,rho] = cart2pol(res(:,1),res(:,2));
theta = theta + pi;
theta = theta * 180/pi;%Se transforma a grados
[sorTheta,pos] = sort(theta);
%rho = rho(pos);

%PASO 2
%LOESS Regression
Ysmooth = malowess(theta,rho,'Order',2);
% y = fLOESS([sorTheta,rho],.2);
figure,plot(sorTheta,Ysmooth)
hold on
plot(sorTheta,rho)

%Se Invierte la forma polar a cartesiana, para obtener la figura, pero con
%suavidad
theta = theta * pi/180;
theta = theta - pi;
[nX,nY] = pol2cart(theta,Ysmooth);
nX = nX + x;
nY = nY + y;
figure,plot(nX,nY);

%PASO 3
%Obtener los puntos del contorno que tengan una distancia euclidiana que no
%rebase un valor epsilon
ep = .3;
parPuntos = [];
first_point = [];
count_par_points = 1;
evaluado = 0;
for j = 1:size(nX,1)
    if(isempty(first_point))
        first_point = [nX(j),nY(j)];
        last_points = first_point;
        count_points = 2;
    else
        second_point = [nX(j),nY(j)];
        last_points(count_points,:) = second_point;
        ec_dist = get_ec_dist(first_point,second_point);
        ac_dist = get_ac_dist(last_points);
        if abs(ec_dist-ac_dist) <= ep
            evaluado = 1;
            count_points = count_points + 1;
        else
            if(evaluado)
                parPuntos(count_par_points,:) = first_point;
                parPuntos(count_par_points+1,:) = second_point;
                count_par_points = count_par_points + 2;
                evaluado = 0;
            end
            first_point = [];
        end
    end
end


puntosEliminar = [];
count = 1;
for  i = 1:2:size(parPuntos,1)
    puntoMedio = [sum(parPuntos(i:i+1,1))/2,sum(parPuntos(i:i+1,2))/2];
    [xi,yi] = polyxpoly([puntoMedio(1),centroid(2)],[puntoMedio(2),centroid(1)],nX,nY); %Se revisa si intersecta la linea entre el punto medio de los 2 puntos y el centroide de la figura
    %Si estas variables traen algo?, ya no cumple porque estan
    %intersectando antes la figura que la linea generada, y se sacan del
    %arreglo
    if(~isempty(xi))
        puntosEliminar(count) = i;
        puntosEliminar(count+1) = i+1;
        count = count +2;
    end
end
parPuntos(puntosEliminar,:) = []; %Se eliminan los puntos que no cumplen
hold on
scatter(parPuntos(:,1),parPuntos(:,2));
scatter(centroid(2),centroid(1)); %Centroide de la figura

for  i = 1:2:size(parPuntos,1)
    plot(parPuntos(i:i+1,1),parPuntos(i:i+1,2))
    puntoMedio = [sum(parPuntos(i:i+1,1))/2,sum(parPuntos(i:i+1,2))/2];
    plot([puntoMedio(1),centroid(2)],[puntoMedio(2),centroid(1)])
end

%Acumulado de la distancia euclidiana de los puntos
function ac_dist = get_ac_dist(points)
    ac_dist = 0;
    for m = 1:size(points,1) - 1
       ac_dist = ac_dist + get_ec_dist(points(m,:),points(m+1,:));
    end
end

%Distancia euclidiana
function ec_dist = get_ec_dist(point1,point2)
    ec_dist = sqrt((point1(1)-point2(1))^2+(point1(2)-point2(2))^2);
end
