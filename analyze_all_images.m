am=57681;



% %% %% for class one
armon=40;
dist_max=10;
pathi=[pwd,'\class1\'];
f=dir([pathi,'*.jpg']);
cont=length(f)+1; 

pathi=[pwd,'\failsize\'];
f=dir([pathi,'*.jpg']);
cont2=length(f)+1;  

pathi=[pwd,'\all\'];
f=dir([pathi,'*.jpg']);
file={};
for i=1:length(f)
file(i)={f(i).name};
end
path=pathi;

roundness=zeros(1);
for i=1:length(file)
    im = imread([path, file{i}]);
    im = imbinarize(im(:,:,1));
    im=regionprops(im,'Image');
    im=im(1).Image;
    area=sum(im(:));
    area=area/am;
    if area>0.6
        try
        r = forclass1([path, file{i}],armon,dist_max,file{i},'results_db9',0);
        if r<0.16
            roundness(cont)=r;
              movefile([path, file{i}], ['C:\Users\Gama\Documents\Github\Morphology\class1\',num2str(cont),'.jpg'])
            cont=cont+1;
        end
        catch
            fprintf('fail')
        end
    else
         movefile([path, file{i}], ['C:\Users\Gama\Documents\Github\Morphology\fail_size\',num2str(cont2),'.jpg'])
        cont2=cont2+1;
    end
    
end

try
r=open('C:\Users\Gama\Documents\Github\Morphology\class1\roundness.mat');
r=r.roundness;
roundness=[roundness,r];
catch
roundness=roundness/max(roundness);
roundness=(1+roundness-mean(roundness))/10;
end
save('C:\Users\Gama\Documents\Github\Morphology\class1\roundness.mat','roundness');

%% % for class two
armon=25;
dist_max=10;
pathi=[pwd,'\class2\'];
f=dir([pathi,'*.jpg']);
cont=length(f)+1; 

pathi=[pwd,'\failsize\'];
f=dir([pathi,'*.jpg']);
cont2=length(f)+1;  

pathi=[pwd,'\all\'];
f=dir([pathi,'*.jpg']);
file={};
for i=1:length(f)
file(i)={f(i).name};
end
path=pathi;

roundness=zeros(1);
for i=1:length(file)
    im = imread([path, file{i}]);
    im = imbinarize(im(:,:,1));
    im=regionprops(im,'Image');
    im=im(1).Image;
    area=sum(im(:));
    area=area/am;
    if area>0.6
        try
        r = forclass2([path, file{i}],armon,dist_max,file{i},'results_db9',0);
        if r<0.309
            roundness(cont)=r;
              movefile([path, file{i}], ['C:\Users\Gama\Documents\Github\Morphology\class2\',num2str(cont),'.jpg'])
            cont=cont+1;
        end
        catch
            fprintf('fail')
        end
    else
         movefile([path, file{i}], ['C:\Users\Gama\Documents\Github\Morphology\fail_size\',num2str(cont2),'.jpg'])
        cont2=cont2+1;
    end
    
end


try
r=open('C:\Users\Gama\Documents\Github\Morphology\class2\roundness.mat');
r=r.roundness;
roundness=[roundness,r];
catch
roundness=roundness/max(roundness);
roundness=(2+roundness-mean(roundness))/10;
end
save('C:\Users\Gama\Documents\Github\Morphology\class2\roundness.mat','roundness');


%% % for class three
armon=40;
dist_max=10;
pathi=[pwd,'\class3\'];
f=dir([pathi,'*.jpg']);
cont=length(f)+1; 

pathi=[pwd,'\failsize\'];
f=dir([pathi,'*.jpg']);
cont2=length(f)+1;  

pathi=[pwd,'\all\'];
f=dir([pathi,'*.jpg']);
file={};
for i=1:length(f)
file(i)={f(i).name};
end
path=pathi;

roundness=zeros(1);
for i=1:length(file)
    im = imread([path, file{i}]);
    im = imbinarize(im(:,:,1));
    im=regionprops(im,'Image');
    im=im(1).Image;
    area=sum(im(:));
    area=area/am;
    if area>0.6
        try
        r = forclass3([path, file{i}],armon,dist_max,file{i},'results_db9',0);
        if r<0.36
            roundness(cont)=r;
              movefile([path, file{i}], ['C:\Users\Gama\Documents\Github\Morphology\class3\',num2str(cont),'.jpg'])
            cont=cont+1;
        end
        catch
            fprintf('fail')
        end
    else
         movefile([path, file{i}], ['C:\Users\Gama\Documents\Github\Morphology\fail_size\',num2str(cont2),'.jpg'])
        cont2=cont2+1;
    end
    
end

try
r=open('C:\Users\Gama\Documents\Github\Morphology\class3\roundness.mat');
r=r.roundness;
roundness=[r,roundness];
catch
roundness=roundness/max(roundness);
roundness=(3+roundness-mean(roundness))/10;
end
save('C:\Users\Gama\Documents\Github\Morphology\class3\roundness.mat','roundness');


%% % for class four
armon=25;
dist_max=15;
pathi=[pwd,'\class4\'];
f=dir([pathi,'*.jpg']);
cont=length(f)+1; 

pathi=[pwd,'\failsize\'];
f=dir([pathi,'*.jpg']);
cont2=length(f)+1;  

pathi=[pwd,'\all\'];
f=dir([pathi,'*.jpg']);
file={};
for i=1:length(f)
file(i)={f(i).name};
end
path=pathi;

roundness=zeros(1);
for i=1:length(file)
    im = imread([path, file{i}]);
    im = imbinarize(im(:,:,1));
    im=regionprops(im,'Image');
    im=im(1).Image;
    area=sum(im(:));
    area=area/am;
    if area>0.6
        try
        r = forclass4([path, file{i}],armon,dist_max,file{i},'results_db9',0);
        if r<0.496
            roundness(cont)=r;
              movefile([path, file{i}], ['C:\Users\Gama\Documents\Github\Morphology\class4\',num2str(cont),'.jpg'])
            cont=cont+1;
        end
        catch
            fprintf('fail')
        end
    else
         movefile([path, file{i}], ['C:\Users\Gama\Documents\Github\Morphology\fail_size\',num2str(cont2),'.jpg'])
        cont2=cont2+1;
    end
    
end
try
r=open('C:\Users\Gama\Documents\Github\Morphology\class4\roundness.mat');
r=r.roundness;
roundness=[r,roundness];
catch
roundness=roundness/max(roundness);
roundness=(4+roundness-mean(roundness))/10;
end
save('C:\Users\Gama\Documents\Github\Morphology\class4\roundness.mat','roundness');


%% % for class five
armon=25;
dist_max=15;
pathi=[pwd,'\class5\'];
f=dir([pathi,'*.jpg']);
cont=length(f)+1; 

pathi=[pwd,'\failsize\'];
f=dir([pathi,'*.jpg']);
cont2=length(f)+1;  

pathi=[pwd,'\all\'];
f=dir([pathi,'*.jpg']);
file={};
for i=1:length(f)
file(i)={f(i).name};
end
path=pathi;

roundness=zeros(1);
for i=1:length(file)
    im = imread([path, file{i}]);
    im = imbinarize(im(:,:,1));
    im=regionprops(im,'Image');
    im=im(1).Image;
    area=sum(im(:));
    area=area/am;
    if area>0.6
        try
        r = forclass5([path, file{i}],armon,dist_max,file{i},'results_db9',0);
        if r<0.61
            roundness(cont)=r;
              movefile([path, file{i}], ['C:\Users\Gama\Documents\Github\Morphology\class5\',num2str(cont),'.jpg'])
            cont=cont+1;
        end
        catch
            fprintf('fail')
        end
    else
         movefile([path, file{i}], ['C:\Users\Gama\Documents\Github\Morphology\fail_size\',num2str(cont2),'.jpg'])
        cont2=cont2+1;
    end
    
end


try
r=open('C:\Users\Gama\Documents\Github\Morphology\class5\roundness.mat');
r=r.roundness;
roundness=[r,roundness];
catch
roundness=roundness/max(roundness);
roundness=(5+roundness-mean(roundness))/10;
end
save('C:\Users\Gama\Documents\Github\Morphology\class5\roundness.mat','roundness');



%% % for class SIX
armon=35;
dist_max=10;
pathi=[pwd,'\class6\'];
f=dir([pathi,'*.jpg']);
cont=length(f)+1; 

pathi=[pwd,'\failsize\'];
f=dir([pathi,'*.jpg']);
cont2=length(f)+1;  

pathi=[pwd,'\all\'];
f=dir([pathi,'*.jpg']);
file={};
for i=1:length(f)
file(i)={f(i).name};
end
path=pathi;

roundness=zeros(1);
for i=1:length(file)
    im = imread([path, file{i}]);
    im = imbinarize(im(:,:,1));
    im=regionprops(im,'Image');
    im=im(1).Image;
    area=sum(im(:));
    area=area/am;
    if area>0.6
        try
        r = forclass6([path, file{i}],armon,dist_max,file{i},'results_db9',0);
        if r<0.639
            roundness(cont)=r;
              movefile([path, file{i}], ['C:\Users\Gama\Documents\Github\Morphology\class6\',num2str(cont),'.jpg'])
            cont=cont+1;
        end
        catch
            fprintf('fail')
        end
    else
         movefile([path, file{i}], ['C:\Users\Gama\Documents\Github\Morphology\fail_size\',num2str(cont2),'.jpg'])
        cont2=cont2+1;
    end
    
end
try
r=open('C:\Users\Gama\Documents\Github\Morphology\class6\roundness.mat');
r=r.roundness;
roundness=[r,roundness];
catch
roundness=roundness/max(roundness);
roundness=(6+roundness-mean(roundness))/10;
end
save('C:\Users\Gama\Documents\Github\Morphology\class6\roundness.mat','roundness');


%% % for class SEVEN
armon=25;
dist_max=15;
pathi=[pwd,'\class7\'];
f=dir([pathi,'*.jpg']);
cont=length(f)+1; 

pathi=[pwd,'\failsize\'];
f=dir([pathi,'*.jpg']);
cont2=length(f)+1;  

pathi=[pwd,'\all\'];
f=dir([pathi,'*.jpg']);
file={};
for i=1:length(f)
file(i)={f(i).name};
end
path=pathi;

roundness=zeros(1);
for i=1:length(file)
    im = imread([path, file{i}]);
    im = imbinarize(im(:,:,1));
    im=regionprops(im,'Image');
    im=im(1).Image;
    area=sum(im(:));
    area=area/am;
    if area>0.6
        try
        r = forclass7([path, file{i}],armon,dist_max,file{i},'results_db9',0);
        if r<0.843
            roundness(cont)=r;
              movefile([path, file{i}], ['C:\Users\Gama\Documents\Github\Morphology\class7\',num2str(cont),'.jpg'])
            cont=cont+1;
        end
        catch
            fprintf('fail')
        end
    else
         movefile([path, file{i}], ['C:\Users\Gama\Documents\Github\Morphology\fail_size\',num2str(cont2),'.jpg'])
        cont2=cont2+1;
    end
    
end
try
r=open('C:\Users\Gama\Documents\Github\Morphology\class7\roundness.mat');
r=r.roundness;
roundness=[r,roundness];
catch
roundness=roundness/max(roundness);
roundness=(7+roundness-mean(roundness))/10;
end
save('C:\Users\Gama\Documents\Github\Morphology\class7\roundness.mat','roundness');


%% % for class eigth
armon=25;
dist_max=15;
pathi=[pwd,'\class8\'];
f=dir([pathi,'*.jpg']);
cont=length(f)+1; 

pathi=[pwd,'\failsize\'];
f=dir([pathi,'*.jpg']);
cont2=length(f)+1;  

pathi=[pwd,'\all\'];
f=dir([pathi,'*.jpg']);
file={};
for i=1:length(f)
file(i)={f(i).name};
end
path=pathi;

roundness=zeros(1);
for i=1:length(file)
    im = imread([path, file{i}]);
    im = imbinarize(im(:,:,1));
    im=regionprops(im,'Image');
    im=im(1).Image;
    area=sum(im(:));
    area=area/am;
    if area>0.6
        try
        r = forclass8([path, file{i}],armon,dist_max,file{i},'results_db9',0);
        if r<0.94
            roundness(cont)=r;
              movefile([path, file{i}], ['C:\Users\Gama\Documents\Github\Morphology\class8\',num2str(cont),'.jpg'])
            cont=cont+1;
        end
        catch
            fprintf('fail')
        end
    else
         movefile([path, file{i}], ['C:\Users\Gama\Documents\Github\Morphology\fail_size\',num2str(cont2),'.jpg'])
        cont2=cont2+1;
    end
    
end

try
r=open('C:\Users\Gama\Documents\Github\Morphology\class8\roundness.mat');
r=r.roundness;
roundness=[r,roundness];
roundness=roundness/max(roundness);
roundness=8+roundness;
catch
roundness=roundness/max(roundness);
roundness=(8+roundness-mean(roundness))/10;
end
save('C:\Users\Gama\Documents\Github\Morphology\class8\roundness.mat','roundness');





% armon=80;
% dist_max=10;
% 
% 
% pathi=[pwd,'\class1\'];
% f=dir([pathi,'*.jpg']);
% file={};
% for i=1:length(f)
%     file(i)={f(i).name};
% end
% path=pathi;
% 
% 
% roundness=zeros(1);
% for i=1:length(file)
%     im = imread([path, file{i}]);
%     im = imbinarize(im(:,:,1));
%     im=regionprops(im,'Image');
%     im=im(1).Image;
%     roundness(i) = inscribedCircles_6 ([path, file{i}],armon,dist_max,file{i},'results_db9',0);
%     movefile([path, file{i}], ['C:\Users\Gama\Documents\Github\Morphology\class1t\',num2str(i),'.jpg'])    
%     
% end
% 



%%%%%for class one 




% if size(results,1)>6
% % [~,ind]= min([results{:,1}]);
% % results(ind,:)=[];
% [~,ind]= max([results{:,1}]);
% results(ind,:)=[];
% [~,ind]= max([results{:,1}]);
% results(ind,:)=[];
% % [~,ind]= min([results{:,1}]);
% % results(ind,:)=[];
% [~,ind]= max([results{:,1}]);
% results(ind,:)=[];
% end


%if gof.rsquare<0.90
% R<=radii










%%%%%for class two

% 
% [results]=fitting_circ(xt,yt,radii,dist_max,imSmooth,0);
% if size(results,1)>6
%     [~,oo]=sort([results{:,1}],'descend');
%     results=results(oo(1:6),:);
% % [~,ind]= min([results{:,1}]);
% % % results(ind,:)=[];
% % [~,ind]= max([results{:,1}]);
% % results(ind,:)=[];
% % [~,ind]= max([results{:,1}]);
% % results(ind,:)=[];
% % % [~,ind]= min([results{:,1}]);
% % % results(ind,:)=[];
% % [~,ind]= max([results{:,1}]);
% % results(ind,:)=[];
% % [~,ind]= max([results{:,1}]);
% % results(ind,:)=[];
% % [~,ind]= max([results{:,1}]);
% % results(ind,:)=[];
% end
% 







