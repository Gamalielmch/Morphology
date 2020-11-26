%Lectura de la imagen
I=imread('C:\Users\Gama\Documents\Github\Morphology\images_test\a1.jpg');

%%Binaraizaci�n de la imagen
bw=im2bw(I(:,:,1)); 

%rotaci�n para visualizaci�n 
bw=fliplr(imrotate(bw,90)); 

% N�mero de armonicos de la transformada de Fourier El�ptico
nc=40;  

%Transformada de Fourier El�ptico
[a,b,c,d,A0,C0,mn,dx,dy,xo,yo]=elifous(bw,nc,0);
% [cons,excen] = armonicosFourierEliptico(nc,'C:\Users\Gama\Documents\Github\Morphology\images_test\a1.jpg');

%magnitud y fase
Mx=sqrt(a.^2+b.^2);
My=sqrt(c.^2+d.^2);
Ax=unwrap(atan(b./a));
Ay=unwrap(atan(d./c));

%Normalizaci�n de escala de la part�cula
Mx=Mx/max(Mx);
My=My/max(My);

%Estandarizaci�n z-score
Mxs= (Mx(2:end)-mean(Mx(2:end)))./std(Mx(2:end));
Mys= (My(2:end)-mean(My(2:end)))./std(My(2:end));

%Renormalizando 
Mxs=Mxs-min(Mxs);
Mxs=Mxs/max(Mxs);

Mys=Mys-min(Mys);
Mys=Mys/max(Mys);



