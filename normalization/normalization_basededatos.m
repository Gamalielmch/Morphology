%Lectura de la imagen
dataset=load('db_roudness.mat');
spectra=dataset.spectra;
roundness=dataset.roundness;
Mx=zeros(size(spectra,1),39);
My=zeros(size(spectra,1),39);

for i=1:size(spectra,1)
a=spectra(i,1:40);
b=spectra(i,41:80);
c=spectra(i,81:120);
d=spectra(i,121:160);
Mxt=sqrt(a.^2+b.^2);
Myt=sqrt(c.^2+d.^2);

%Normalización de escala de la partícula
Mxt=Mxt/max(Mxt);
Myt=Myt/max(Myt);
%Estandarización z-score
Mxs= (Mxt(2:end)-mean(Mxt(2:end)))./std(Mxt(2:end));
Mys= (Myt(2:end)-mean(Myt(2:end)))./std(Myt(2:end));

%Renormalizando 
Mxs=Mxs-min(Mxs);
Mx(i,:)=Mxs/max(Mxs);

Mys=Mys-min(Mys);
My(i,:)=Mys/max(Mys);
end








