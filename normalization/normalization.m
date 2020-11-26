%Lectura de la imagen
I=imread('C:\Users\Gama\Documents\Github\Morphology\images_test\a1.jpg');

%%Binaraización de la imagen
bw=im2bw(I(:,:,1)); 

%rotación para visualización 
bw=fliplr(imrotate(bw,90)); 

% Número de armonicos de la transformada de Fourier Elíptico
nc=40;  

%Transformada de Fourier Elíptico
[a,b,c,d,A0,C0,mn,dx,dy,xo,yo]=elifous(bw,nc,0);
% [cons,excen] = armonicosFourierEliptico(nc,'C:\Users\Gama\Documents\Github\Morphology\images_test\a1.jpg');

%magnitud y fase
Mx=sqrt(a.^2+b.^2);
My=sqrt(c.^2+d.^2);
Ax=unwrap(atan(b./a));
Ay=unwrap(atan(d./c));

%Normalización de escala de la partícula
Mx=Mx/max(Mx);
My=My/max(My);

%Estandarización z-score
Mxs= (Mx(2:end)-mean(Mx(2:end)))./std(Mx(2:end));
Mys= (My(2:end)-mean(My(2:end)))./std(My(2:end));

%Renormalizando 
Mxs=Mxs-min(Mxs);
Mxs=Mxs/max(Mxs);

Mys=Mys-min(Mys);
Mys=Mys/max(Mys);



