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

roundness=zeros(9,9);
for j=1:9
    for i=1:9
        roundness(j,i) = inscribedCircles_5 ([path, file{i+(j-1)*9}],armon,dist_max,file{i+(j-1)*9},0);
    end
end
color=lines(9);
figure
plot(1:9,roundness(1,:),'s','color',color(1,:),'MarkerFaceColor',color(1,:))
hold on
plot(10:18,roundness(2,:),'s','color',color(2,:),'MarkerFaceColor',color(2,:))
plot(19:27,roundness(3,:),'s','color',color(3,:),'MarkerFaceColor',color(3,:))
plot(28:36,roundness(4,:),'s','color',color(4,:),'MarkerFaceColor',color(4,:))
plot(37:45,roundness(5,:),'s','color',color(5,:),'MarkerFaceColor',color(5,:))
plot(46:54,roundness(6,:),'s','color',color(6,:),'MarkerFaceColor',color(6,:))
plot(55:63,roundness(7,:),'s','color',color(7,:),'MarkerFaceColor',color(7,:))
plot(64:72,roundness(8,:),'s','color',color(8,:),'MarkerFaceColor',color(8,:))
plot(73:81,roundness(9,:),'s','color',color(9,:),'MarkerFaceColor',color(9,:))


