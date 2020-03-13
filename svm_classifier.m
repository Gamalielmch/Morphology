%%% Frecuencias bunas: 0.0166 0.0431


%class well rounded fc=0.024
%class rounded fc=0.04691
%class sub-rounded fc=0.1404
%class angular fc=0.1404

clear
padi='C:\Users\Usuario\Documents\MATLAB\particulas';
[filenamed, pathd]= uigetfile({'*.*'},'Images of CLASS ','MultiSelect', 'on',padi);
if pathd==0
    error('No files, try again, run function again');
end

assi=[ones(1,5), ones(1,5)*2, ones(1,5)*3, ones(1,5)*4, ones(1,5)*5]';

if iscell(filenamed)
    nclass=length(filenamed);
    images=cell(nclass,1);
    for i=1:nclass
        I=imread([pathd, filenamed{i}]);
        if size(I,3)>3
            I=I(:,:,1:3);
        end
        images{i}=im2bw(mat2gray(I));
        
    end
else
    I=imread([pathd, filenamed]);
    if size(I,3)>3
        I=I(:,:,1:3);
    end
    images{1}=im2bw(mat2gray(I));
    nclass=1;
end
nc=128;   %%% número de armónicos
ncr=3;    %%% Número de armónicos para la reconstrucción
graf=1;   %%% visualización de la reconstrucción
ii = 1: nc;
% figure, hold on
a=zeros(nclass,nc);
b=zeros(nclass,nc);
c=zeros(nclass,nc);
d=zeros(nclass,nc);
mn=zeros(nc,1);
for i=1:nclass
    
    [a(i,:),b(i,:),c(i,:),d(i,:),~,~,mn(i)]=elifou(images{i},nc,1);
    realt=zeros(mn(i),2);
    for t = 1 : mn(i)
        realt(t,1) =sum( (a(i,ii)  .* cos(2 .* ii .* pi .* t / mn(i)) + b(i,ii)  .* sin(2 .* ii .* pi .* t ./ mn(i))));
        realt(t,2) = sum((c(i,ii)  .* cos(2 .* ii .* pi .* t / mn(i)) + d(i,ii)  .* sin(2 .* ii .* pi .* t ./ mn(i))));
    end
    %     plot(realt(:,1), realt(:,2), 'b', 'linewidth', 3);
    real{i}={realt};
end



fc1=0.0166;
fc2=0.0431;

Response1=zeros(1, 25);
Response2=zeros(1, 25);


coef= fir1(200,fc1,'low');
[h,~] = freqz(coef,1,length(a));
h1=abs(h)';

coef= fir1(200,fc2,'low');
[h,~] = freqz(coef,1,length(a));
h2=abs(h)';

for i=1:nclass
    recon=zeros(mn(i),2);
    recon2=zeros(mn(i),2);
    recon(:,1)=sum(h1.*(a(i,ii).* cos(2 *  ii.* pi .* [1 : mn(i)]' / mn(i)) + b(i,ii).* sin(2 * ii.* pi .* [1 : mn(i)]' / mn(i)) ),2  );
    recon(:,2)=sum(h1.*(c(i,ii).* cos(2 *  ii.* pi .* [1 : mn(i)]' / mn(i)) + d(i,ii).* sin(2 * ii.* pi .* [1 : mn(i)]' / mn(i)) ),2  );
    recon2(:,1)=sum(h2.*(a(i,ii).* cos(2 *  ii.* pi .* [1 : mn(i)]' / mn(i)) + b(i,ii).* sin(2 * ii.* pi .* [1 : mn(i)]' / mn(i)) ),2  );
    recon2(:,2)=sum(h2.*(c(i,ii).* cos(2 *  ii.* pi .* [1 : mn(i)]' / mn(i)) + d(i,ii).* sin(2 * ii.* pi .* [1 : mn(i)]' / mn(i)) ),2  );
    realt=real{i};
    realt=realt{1};
    Response1(i)=sum(   sqrt(  (recon(:,1)- realt(:,1)).^2 +  (recon(:,2)- realt(:,2)).^2  )  ) /...
        sum(   sqrt(  (recon(:,1)).^2 +  (recon(:,2)).^2  )  );
    Response2(i)=sum(   sqrt(  (recon2(:,1)- realt(:,1)).^2 +  (recon2(:,2)- realt(:,2)).^2  )  ) /...
        sum(   sqrt(  (recon2(:,1)).^2 +  (recon2(:,2)).^2  )  );
end

figure
plot(Response1(1:5),Response2(1:5),'r+')
hold on
plot(Response1(6:10),Response2(6:10),'g+')
plot(Response1(11:15),Response2(11:15),'b+')
plot(Response1(16:20),Response2(16:20),'k+')
plot(Response1(21:25),Response2(21:25),'c+')


label=cell(25,1);
for i=1:25
    if i<6
label(i)={'Angular'};
    elseif i>5 && i<11
      label(i)={'Sub-angular'};
    elseif i>10 && i<16
      label(i)={'Sub-rounded'};
    elseif i>15 && i<21 
     label(i)={'Rounded'};
    else
    label(i)={'well-rounded'};
    end 
end


figure
gscatter(Response1,Response2,label);
h = gca;
lims = [h.XLim h.YLim]; % Extract the x and y axis limits
title('{\bf Scatter Diagram of Response to filters}');
xlabel('R^2 Filter 1');
ylabel('R^2 Filter 2');
legend('Location','Northwest');

classes = unique(label);
SVMModels = cell(5,1);

for j = 1:numel(classes)
indx = strcmp(label,classes(j)); % Create binary classes for each classifier
SVMModels{j} = fitcsvm([Response1',Response2'],indx,'ClassNames',[false true],'Standardize',true,...
'KernelFunction','rbf','BoxConstraint',1);
end


d = 0.0001;
[x1Grid,x2Grid] = meshgrid(min(Response1):d:max(Response1),...
    min(Response2):d:max(Response2));
xGrid = [x1Grid(:),x2Grid(:)];
N = size(xGrid,1);




Scores = zeros(N,numel(classes));
for j = 1:numel(classes)
    [~,score] = predict(SVMModels{j},xGrid);
    Scores(:,j) = score(:,2); % Second column contains positive-class scores
end

[~,maxScore] = max(Scores,[],2);

col= lines(5);
figure
h(1:5) = gscatter(xGrid(:,1),xGrid(:,2),maxScore,col);
hold on
h(6:10) = gscatter(Response1,Response2,label);
title('{\bf Roundness Classification Regions}');
xlabel('R^2 Filter 1');
ylabel('R^2 Filter 2');
legend(h,{'Angular region','Sub-angular region','Sub-rounded region','Rounded region','Well-rounded region',...
    'Observed angular','Observed Sub-angular','Observed Sub-rounded','Observed Rounded ','Observed Well-rounded'},...
    'Location','Northwest');
axis tight
hold off




