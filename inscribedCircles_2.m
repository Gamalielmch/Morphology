function main
close all;
%load franke
im = imread('testimage/e1.jpg');
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

% res = boundary - [x,y]; %Se resta las posiciones del centro del círuclo mas grande, para poner la figura en el origen
% [theta,rho] = cart2pol(res(:,1),res(:,2));
% theta = theta + pi;
% theta = theta * 180/pi;%Se transforma a grados
% [sorTheta,pos] = sort(theta);
%rho = rho(pos);

%PASO 2
%LOESS Regression
%Ysmooth = smooth(theta,rho,0.2,'rloess');
imSmooth = bw_chain_fil(BW,30,0); % Fourier eliptico smooth
%imSmooth(size(imSmooth,1)+1,:) = imSmooth(1,:); %Para eliminar el hueco que se genera entre el primero y el ultimo punto, solo se agrega el primer punto al final
% Ysmooth = malowess(theta,rho,'Order',2);
% Ysmooth = smooth(theta,Ysmooth,0.2,'moving');
% % y = fLOESS([sorTheta,rho],.2);
% figure,plot(sorTheta,Ysmooth)
% hold on
plot(imSmooth(:,2),imSmooth(:,1))

%Se Invierte la forma polar a cartesiana, para obtener la figura, pero con
%suavidad
% theta = theta * pi/180;
% theta = theta - pi;
% [nX,nY] = pol2cart(theta,Ysmooth);
% [nXr,nYr] = pol2cart(theta,rho);
% nX = nX + x;
% nY = nY + y;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PASO 3
%Obtener puntos que rebasen un delta
len = size(imSmooth,1); %numero de pixeles en el contorno
puntos = 1;
i = 1;
cont = 1;
delta = 0.4;
%figure, plot(Ysmooth)
%hold on
%plot(i,Ysmooth(i),'ks')
while i<len
    for xi = len:-1:i
        c = polyfit([i,xi],[imSmooth(i,2),imSmooth(xi,2)],1);
        ylin = polyval(c,i:xi)';
        maxi_div=max(abs(imSmooth(i:xi,2)-ylin));
        if maxi_div<delta
            %plot(xi,imSmooth(xi,2),'ks')
            %plot([i,xi],[ylin(1),ylin(end)])
            %
            cont=cont+1;
            puntos(cont)=xi;
            i = xi;
            break
        end
    end
end
%figx=figure;
plot(imSmooth(:,2),imSmooth(:,1))
hold on
for i=1:length(puntos)
   plot(imSmooth(puntos(i),2),imSmooth(puntos(i),1),'ks') 
   drawnow
end
for  i = 1:length(puntos)-1
    plot([imSmooth(puntos(i),2),imSmooth(puntos(i+1),2)],[imSmooth(puntos(i),1),imSmooth(puntos(i+1),1)])
    %puntoMedio = [sum(parPuntos(i:i+1,1))/2,sum(parPuntos(i:i+1,2))/2];
    %plot([puntoMedio(1),centroid(2)],[puntoMedio(2),centroid(1)])
end
puntosSeleccionados = [];
count = 1;

for  i = 1:length(puntos)-1
    puntoMedio = [sum(imSmooth(puntos(i:i+1),1))/2,sum(imSmooth(puntos(i:i+1),2))/2];
    [xi,yi] = polyxpoly(imSmooth(puntos(i:i+1),1),imSmooth(puntos(i:i+1),2),imSmooth(:,1),imSmooth(:,2)); %Se revisa si intersecta la linea entre el punto medio de los 2 puntos y el centroide de la figura
    [xi2,~] = polyxpoly([puntoMedio(1),centroid(1)],[puntoMedio(2),centroid(2)],imSmooth(:,1),imSmooth(:,2));
    %plot(yi,xi,'md');
    %Si estas variables traen algo?, ya no cumple porque estan
    %intersectando antes la figura que la linea generada, y se sacan del
    xi([1,length(xi)]) = [];
    %arreglo
    if(isempty(xi) && isempty(xi2))
        %if(puntosEliminar.ismember()
        if(~ any(puntosSeleccionados == i)) 
            puntosSeleccionados(count) = i;
            count = count + 1;
        end
        
        if(~ any(puntosSeleccionados == i+1)) 
            puntosSeleccionados(count) = i + 1;
            count = count + 1;
        end
    end
end
figure,plot(imSmooth(:,2),imSmooth(:,1))
puntos = puntos(puntosSeleccionados); %Se eliminan los puntos que no cumplen
hold on
plot(imSmooth(puntos,2),imSmooth(puntos,1),'ks');
title("Puntos seleccionados")
%scatter(centroid(2),centroid(1)); %Centroide de la figura

% for  i = 1:2:length(puntos)
%     plot([imSmooth(puntos(i),2),imSmooth(puntos(i+1),2)],[imSmooth(puntos(i),1),imSmooth(puntos(i+1),1)])
%     %puntoMedio = [sum(parPuntos(i:i+1,1))/2,sum(parPuntos(i:i+1,2))/2];
%     %plot([puntoMedio(1),centroid(2)],[puntoMedio(2),centroid(1)])
% end
end% main

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
