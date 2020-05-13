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
armon=35;
dist_max=15;

roundness=zeros(5,5);
for j=1:5
    for i=1:5
        roundness(j,i) = inscribedCircles_5 ([path, file{i+(j-1)*5}],armon,dist_max,file{i+(j-1)*5},1);
    end
end
color=lines(5);
figure
plot(1:5,roundness(1,:),'s','color',color(1,:),'MarkerFaceColor',color(1,:))
hold on
plot(6:10,roundness(2,:),'s','color',color(2,:),'MarkerFaceColor',color(2,:))
plot(11:15,roundness(3,:),'s','color',color(3,:),'MarkerFaceColor',color(3,:))
plot(16:20,roundness(4,:),'s','color',color(4,:),'MarkerFaceColor',color(4,:))
plot(21:25,roundness(5,:),'s','color',color(5,:),'MarkerFaceColor',color(5,:))


