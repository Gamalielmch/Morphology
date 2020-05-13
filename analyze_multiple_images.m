% function analyze_multiple_images
%%%Read multiple image

pathi=[pwd,'\database1'];
[file,path] = uigetfile({'*.jpg;*.png;*.bmp;*.jpeg;'},'Select One or More Files','MultiSelect', 'on',pathi);
if iscell(file)
    ni=length(file);
else
    ni=1;
    file={file};
end
armon=35;
dist_max=15;


roundness=zeros(1,ni);
for i=1:ni
roundness(i) = inscribedCircles_5 ([path, file{i}],armon,dist_max,file{i},1);
end


% figure
% plot(1:length(roundness),roundness,'s','color',[1 0 0],'MarkerFaceColor',[1,0,0])
% 
% %%%pettijohn
% class1=roundness(1:5);
% class2=roundness(6:10);
% class3=roundness(11:15);
% class4=roundness(16:20);
% class5=roundness(21:25);
% figure
% plot(1:5,class1,'s','color',[1 0 0],'MarkerFaceColor',[1,0,0])
% hold on
% plot(6:10,class2,'s','color',[0 1 0],'MarkerFaceColor',[0,1,0])
% plot(11:15,class3,'s','color',[0 0 1],'MarkerFaceColor',[0,0,1])
% plot(16:20,class4,'s','color',[1 0 1],'MarkerFaceColor',[1,0,1])
% plot(21:25,class5,'s','color',[0 1 1],'MarkerFaceColor',[0,1,1])
% 
% %%%otras
% figure
% plot(1:10,roundness1,'s','color',[1 0 0],'MarkerFaceColor',[1,0,0])
% hold on
% plot(11:20,roundness2,'s','color',[0 1 0],'MarkerFaceColor',[0,1,0])
% plot(21:30,roundness3,'s','color',[0 0 1],'MarkerFaceColor',[0,0,1])
% plot(31:40,roundness4,'s','color',[1 0 1],'MarkerFaceColor',[1,0,1])
% plot(41:50,roundness5,'s','color',[0 1 1],'MarkerFaceColor',[0,1,1])


% %%%Krumbein
% 
% color=lines(9);
% figure
% plot(1:9,roundness1,'s','color',color(1,:),'MarkerFaceColor',color(1,:))
% hold on
% plot(10:18,roundness2,'s','color',color(2,:),'MarkerFaceColor',color(2,:))
% plot(19:27,roundness3,'s','color',color(3,:),'MarkerFaceColor',color(3,:))
% plot(28:36,roundness4,'s','color',color(4,:),'MarkerFaceColor',color(4,:))
% plot(37:45,roundness5,'s','color',color(5,:),'MarkerFaceColor',color(5,:))
% plot(46:54,roundness6,'s','color',color(6,:),'MarkerFaceColor',color(6,:))
% plot(55:63,roundness7,'s','color',color(7,:),'MarkerFaceColor',color(7,:))
% plot(64:72,roundness8,'s','color',color(8,:),'MarkerFaceColor',color(8,:))
% plot(73:81,roundness9,'s','color',color(9,:),'MarkerFaceColor',color(9,:))


