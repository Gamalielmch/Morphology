% function analyze_multiple_images
%%%Read multiple image

pathi=[pwd,'\testimage'];
[file,path] = uigetfile({'*.jpg;*.png;*.bmp;*.jpeg;'},'Select One or More Files','MultiSelect', 'on',pathi);
if iscell(file)
    ni=length(file);
else
    ni=1;
    file={file};
end
armon=40;
dist_max=15;

%Bueno 2
% armon=40;
% dist_max=15;

%Bueno 1
% armon=30;
% dist_max=15;
roundness=zeros(1,ni);
for i=1:ni
roundness(i) = inscribedCircles_5 ([path, file{i}],armon,dist_max,file{i});
end
class1=roundness(1:5);
class2=roundness(6:10);
class3=roundness(11:15);
class4=roundness(16:20);
class5=roundness(21:25);
figure
plot(1:5,class1,'s','color',[1 0 0],'MarkerFaceColor',[1,0,0])
hold on
plot(6:10,class2,'s','color',[0 1 0],'MarkerFaceColor',[0,1,0])
plot(11:15,class3,'s','color',[0 0 1],'MarkerFaceColor',[0,0,1])
plot(16:20,class4,'s','color',[1 0 1],'MarkerFaceColor',[1,0,1])
plot(21:25,class5,'s','color',[0 1 1],'MarkerFaceColor',[0,1,1])