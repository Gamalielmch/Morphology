% 
% 
% cont=1;
% for jj=1:14
%     
%     pathi=[pwd,'\database',num2str(jj),'\'];
%     f=dir([pathi,'*.jpg']);
%     file={};
%     for i=1:length(f)
%         file(i)={f(i).name};
%     end
%     path=pathi;
%     
%    
%     for i=1:length(file)
%         im = imread([path, file{i}]);
%         im = imbinarize(im(:,:,1));
%         im=regionprops(im,'Image');
%         im=im(1).Image;
%         imwrite(padarray(im,[20 20],'both'),['C:\Users\Gama\Documents\Github\Morphology\all\',num2str(cont),'.jpg']);
%         cont=cont+1;
%     end
% end





% 
% am=57681;
% pathi=[pwd,'\all\'];
% f=dir([pathi,'*.jpg']);
% file={};
% for i=1:length(f)
% file(i)={f(i).name};
% end
% path=pathi;
% 
% for i=1:length(file)
%     im = imread([path, file{i}]);
%     im = imbinarize(im(:,:,1));
%     im=regionprops(im,'Image');
%     im=im(1).Image;
%     im=padarray(im,[20 20],'both');
%     area=sum(im(:));
%     area=area/am;
%     if area<0.61
%     im=imresize(im,sqrt(1/area));
%     im=regionprops(im,'Image');
%     im=im(1).Image;
%     im=padarray(im,[20 20],'both');
%     imwrite(im, ['C:\Users\Gama\Documents\Github\Morphology\all\', file{i}]);
%     end
% %      movefile([path, file{i}], ['C:\Users\Gama\Documents\Github\Morphology\all\', file{i}])
%      
% end



cont=1;
for jj=1:9
    
    pathi=[pwd,'\class',num2str(jj),'\'];
    f=dir([pathi,'*.jpg']);
    file={};
    for i=1:length(f)
        file(i)={f(i).name};
    end
    path=pathi;
    
   
    for i=1:length(file)
        im = imread([path, file{i}]);
%         im = imbinarize(im(:,:,1));
%         im=regionprops(im,'Image');
%         im=im(1).Image;
        imwrite(im,['C:\Users\Gama\Documents\Github\Morphology\allclass\',num2str(cont),'.jpg']);
        cont=cont+1;
    end
end