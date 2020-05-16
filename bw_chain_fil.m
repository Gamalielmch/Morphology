function recon=bw_chain_fil(BW,nc,graf)
%%función basada en Fourier elíptico desarrollado por Kuhl & Giardina
%
%          varibles de entrada:
%  BW: imagen binaria del contorno a analizar, 1 para contorno y 0 para
%  fondo, la imagen no debe contener "holes", debe ser un contorno cerrado
%  nc: Número de coeficientes de las series de Fourier
%
%          varibles de salida:
%  recon= perimetro filtrado


%% Convierte la imagen a tipo binario y realiza ajuste para código de cadena
if ~islogical(BW)
    if size(BW,3)>3
        BW=BW(:,:,1:3);
    end
    BW=im2bw(BW);
end

BW = padarray(BW,2);
BW=bwperim(BW);
BW=bwmorph(BW,'spur');
BW=bwmorph(BW,'thin');
pxl= regionprops(BW,'PixelList');
pxl=pxl.PixelList;
if graf
figure, plot( pxl(:,1),pxl(:,2),'.'), hold on
end
%%                               Código de cadena
[yi,xi]=find(BW,1,'first');
xo=xi; yo=yi;
con=2;
longitud=length(pxl);
BW=double(BW);
le=sqrt(2);

cadena=zeros(1,longitud);
for i=1:longitud-1
    
    pos=[yi,xi+1;yi-1,xi+1;yi-1,xi;yi-1,xi-1;yi,xi-1;yi+1,xi-1;yi+1,xi;yi+1,xi+1];
    if((all(pos(:,1)>0) && all(pos(:,2)>0)) && (all(pos(:,1)<=size(BW,1)) && all(pos(:,2)<=size(BW,2))))
        
        neigh=[BW(pos(1,1),pos(1,2)),BW(pos(2,1),pos(2,2))*le,BW(pos(3,1),pos(3,2)),...
            BW(pos(4,1),pos(4,2))*le,BW(pos(5,1),pos(5,2)),BW(pos(6,1),pos(6,2))*le,...
            BW(pos(7,1),pos(7,2)),BW(pos(8,1),pos(8,2))*le];
    else
        neigh = zeros(1,size(pos,1));
        for m = 1:size(pos,1)
            if(all(pos(m,:)~=0) && all(pos(m,:)<=size(BW)))
               if mod(m,2) == 1
                   neigh(m) = BW(pos(m,1),pos(m,2));
               else
                   neigh(m) = BW(pos(m,1),pos(m,2))*le;
               end
            else
                neigh(m) = 0;
            end
        end
    end
    neigh(neigh==0)=2;
    [~,mino]=min(neigh);
    cadena(i)= mino(1)-1;
    BW(yi,xi)=0;
    xi=pos(mino,2);
    yi=pos(mino,1);
end

pos=[yi,xi+1;yi-1,xi+1;yi-1,xi;yi-1,xi-1;yi,xi-1;yi+1,xi-1;yi+1,xi;yi+1,xi+1];
BW(yo,xo)=1;
neigh=[BW(pos(1,1),pos(1,2)),BW(pos(2,1),pos(2,2))*le,BW(pos(3,1),pos(3,2)),...
    BW(pos(4,1),pos(4,2))*le,BW(pos(5,1),pos(5,2)),BW(pos(6,1),pos(6,2))*le,...
    BW(pos(7,1),pos(7,2)),BW(pos(8,1),pos(8,2))*le];
neigh(neigh==0)=2;
[~,mino]=min(neigh);
cadena(i+1)= mino(1)-1;
mn=length(cadena);


%%                   Tiempo y distancia Traversal
deltat=zeros(1,length(cadena));
deltat(1)=1+((sqrt(2)-1)/2)*(1-(-1)^cadena(1));
deltat2=1+((sqrt(2)-1)/2)*(1-(-1).^cadena);
deltaxa=sign(6-cadena(1))*sign(2-cadena(1));
deltaya=sign(4-cadena(1))*sign(cadena(1));
for i=2:length(cadena)
    deltat(i)=deltat(i-1)+1+((sqrt(2)-1)/2)*(1-(-1)^cadena(i));
    deltaxa(i)=deltaxa(i-1)+ sign(6-cadena(i))*sign(2-cadena(i));
    deltaya(i)=deltaya(i-1)+ sign(4-cadena(i))*sign(cadena(i));
end
T = deltat(end);
deltax=sign(6-cadena).*sign(2-cadena);
deltay=sign(4-cadena).*sign(cadena);




%%                   Cálculo de los coeficientes
a=zeros(1,nc);
b=a;
c=a;
d=a;

q_x = deltax./ deltat2;
q_y = deltay./ deltat2;
for n=1:nc
    two_n_pi = 2 * n * pi;
    sigma_a = 0;
    sigma_b = 0;
    sigma_c = 0;
    sigma_d = 0;
    tp_prev = 0;
    for p=1:length(cadena)
        sigma_a = sigma_a + q_x(p) * (cos(two_n_pi * deltat(p) / T) - cos(two_n_pi * tp_prev / T));
        sigma_b = sigma_b + q_x(p) * (sin(two_n_pi * deltat(p) / T) - sin(two_n_pi * tp_prev / T));
        sigma_c = sigma_c + q_y(p) * (cos(two_n_pi * deltat(p) / T) - cos(two_n_pi * tp_prev / T));
        sigma_d = sigma_d + q_y(p) * (sin(two_n_pi * deltat(p) / T) - sin(two_n_pi * tp_prev / T));
        tp_prev = deltat(p);
    end
    r = T/(2*n^2*pi^2);
    
    a(n) = r * sigma_a;
    b(n) = r * sigma_b;
    c(n)= r * sigma_c;
    d(n)= r * sigma_d;
end

%%                         Cálculo de la componente directa 
A0 = 0;
C0 = 0;
A0 = A0 + deltax(1) / (2 * deltat2(1)) * (deltat(1))^2 ;
C0 = C0 + deltay(1) / (2 * deltat2(1)) * (deltat(1))^2 ;
for p = 2 : length(cadena)
    
    zeta = deltaxa(p - 1) - deltax(p) / deltat2(p) * deltat(p - 1);
    delta = deltaya(p - 1) - deltay(p) / deltat2(p) * deltat(p - 1);
    A0 = A0 + deltax(p) / (2 * deltat2(p)) * ((deltat(p))^2 - (deltat(p - 1))^2) + zeta * (deltat(p) - deltat(p-1));
    C0 = C0 + deltay(p) / (2 * deltat2(p)) * ((deltat(p))^2 - (deltat(p - 1))^2) + delta * (deltat(p) - deltat(p-1));
    
end
A0 = A0 / T;
C0 = C0 / T;




for t = 1 : mn
    x = xo;
    y = yo;
    for ii = 1: nc
        x = x + (a(ii) * cos(2 * ii * pi * t / mn) + b(ii) * sin(2 * ii * pi * t / mn));
        y = y - (c(ii) * cos(2 * ii * pi * t / mn) + d(ii) * sin(2 * ii * pi * t / mn));
    end
    recon(t,1) = A0 + x;
    recon(t,2) = -C0 + y;
    
end
if graf
plot( recon(:,1),recon(:,2))
end