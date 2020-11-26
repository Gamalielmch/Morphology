%Lectura de la imagen
I=imread('C:\Users\Gama\Documents\Github\Morphology\images_test\a1.jpg');

%%Binaraización de la imagen
bw=im2bw(I(:,:,1)); 

%rotación para visualización 
bw=fliplr(imrotate(bw,90)); 

% Número de armonicos de la transformada de Fourier Elíptico
nc=30;  

%Transformada de Fourier Elíptico
[a,b,c,d,A0,C0,mn,dx,dy,xo,yo]=elifous(bw,nc,0);

%magnitud y fase
Mx=sqrt(a.^2+b.^2);
My=sqrt(c.^2+d.^2);
Ax=unwrap(atan(b./a));
Ay=unwrap(atan(d./c));
%mn es el número de pixels en el contorno


%cada coeficiente genera un fasor elipitco, por ejemplo a(1),b(1),d(1) y
%c(1) generan la primera elipse, que es el "promedio" del contorno cerrado
elip=zeros(mn+1,2,nc);
elip2=zeros(mn+1,2,nc);
% A0=0;
% C0=0;

%ciclo para la reconstrucción
for eli=1:nc
        recon=zeros(mn,2);
        recon2=zeros(mn,2);
        for t = 1 : mn
            x = 0;
            y = 0;
            x2 = 0;
            y2 = 0;
            for i = eli : eli
                %serie senos y cosenos
                x =  (a(i) * cos(2 * i * pi * t / mn) + b(i) * sin(2 * i * pi * t / mn)); 
                y = y + (c(i) * cos(2 * i * pi * t / mn) + d(i) * sin(2 * i * pi * t / mn));
                %magnitud y fas
                x2= x2 + Mx(i) * cos(2 * i * pi * t / mn + Ax(i)) ; 
                y2= y2  + My(i) * cos(2 * i * pi * t / mn + Ay(i)) ; 
            end
            recon(t,1) =  x;
            recon(t,2) =  y;
            recon2(t,1) =  x2;
            recon2(t,2) =  y2;
        end
        recon = [recon; recon(1,1) recon(1,2)]; %% n-elipse 
        recon2 = [recon2; recon2(1,1) recon2(1,2)]; %% n-elipse 
    elip(:,:,eli)=recon;
    elip2(:,:,eli)=recon2;
end

%Combinación lineal de las elipse para reconstrucción
rec=zeros(mn,2);
rec2=zeros(mn,2);
for n=1:mn
    rec(n,1)=sum(elip(n,1,1:nc));
    rec(n,2)=sum(elip(n,2,1:nc));
    rec2(n,1)=sum(elip(n,1,1:nc));
    rec2(n,2)=sum(elip(n,2,1:nc));
end

%%% contorno real 
plot( xo+dx, yo+dy,'.','MarkerSize',2,'color',[0,0,1])
hold on
%%% contorno reconstruido
plot( rec(:,1)+A0+xo, rec(:,2)+C0+yo ,'.','MarkerSize',2,'color',[1,0,0])
plot( rec2(:,1)+A0+xo, rec2(:,2)+C0+yo ,'.','MarkerSize',2,'color',[0,1,0])

