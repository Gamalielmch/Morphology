function vari=varmovil(signal,sw,ov)

si=length(signal);
if ~mod(sw,2)
    sw=sw+1;
end
o=1;
mid=((sw-1)/2);
ov=sw-ov;
if ov<=0
    ov=1;
end
for i=mid:ov:si-mid-1
    vari(o)=std(signal(i-mid+1:i+mid+1));
    o=o+1;
end