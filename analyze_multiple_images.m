function analyze_multiple_images
%%%Read multiple image

pathi=[pwd,'\testimage'];
[file,path] = uigetfile({'*.jpg;*.png;*.bmp;*.jpeg;'},'Select One or More Files','MultiSelect', 'on',pathi);
if iscell(file)
    ni=length(file);
else
    ni=1;
    file={file};
end
armon=30;
dist_max=15;
roundness=zeros(1,ni);
for i=1:ni
roundness(i) = inscribedCircles_5 ([path, file{i}],armon,dist_max,file{i});
end