% function analyze_multiple_images
%%%Read multiple image

pathi=[pwd,'\krumbein2'];
[file,path] = uigetfile({'*.jpg;*.png;*.bmp;*.jpeg;'},'Select One or More Files','MultiSelect', 'on',pathi);
if iscell(file)
    ni=length(file);
else
    ni=1;
    file={file};
end
% %%% for class one
%  armon=40;
%  dist_max=10;

%% for class two
% armon=25;
% dist_max=10;

% %%% for class three
% armon=40;
% dist_max=10;

%%% for class four 
% armon=25;
% dist_max=15;

%%% for class five 
% armon=25;
% dist_max=15;

% %%% for class six 
% armon=35;
% dist_max=10;

% %%% for class seven
% armon=25;
% dist_max=15;


% %%% for class eigth 
armon=25;
dist_max=15;
cont=1;
roundness=zeros(9,9);
for j=1:9
    for i=1:9
        roundness(j,i) = forclass7([path, file{i+(j-1)*9}],armon,dist_max,file{i+(j-1)*9},'results_krumbein',0);
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


