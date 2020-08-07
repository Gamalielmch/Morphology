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





% %
% am=63681;
% %pathi=[pwd,'\all\'];
% pathi='C:\Users\Gama\Documents\new_images\';
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
%    % imwrite(im, ['C:\Users\Gama\Documents\Github\Morphology\all\', file{i}]);
%     end
%     imwrite(im, ['C:\Users\Gama\Documents\new_images\temporal\', num2str(i),'.jpg']);
% %      movefile([path, file{i}], ['C:\Users\Gama\Documents\Github\Morphology\all\', file{i}])
%
% end
%


% cont=1;
% for jj=1:9
%
%     pathi=[pwd,'\class',num2str(jj),'\'];
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
% %         im = imbinarize(im(:,:,1));
% %         im=regionprops(im,'Image');
% %         im=im(1).Image;
%         imwrite(im,['C:\Users\Gama\Documents\Github\Morphology\allclass\',num2str(cont),'.jpg']);
%         cont=cont+1;
%     end
% end




% am=63681;
% %pathi=[pwd,'\all\'];
% pathi='C:\Users\Gama\Documents\new_images\';
% f=dir([pathi,'*.jpg']);
% file={};
% for i=1:length(f)
%     file(i)={f(i).name};
% end
% path=pathi;
% ini=0;
% for i=1:length(file)
%     im = imread([path, file{i}]);
%     im = im2bw(im(:,:,1));
%     stats=regionprops(im,'Image');
%     i
%     length(stats)
%     for j=1:length(stats)
%         j
%         im=stats(j).Image;
%         im=padarray(im,[20 20],'both');
%         area=sum(im(:));
%         area=area/am;
%         if area<0.61
%             im=imresize(im,sqrt(1/area));
%             im=regionprops(im,'Image');
%             im=im(1).Image;
%             im=padarray(im,[20 20],'both');
%             % imwrite(im, ['C:\Users\Gama\Documents\Github\Morphology\all\', file{i}]);
%         end
%         imwrite(im, ['C:\Users\Gama\Documents\new_images\database2\', num2str(j+ini),'.jpg']);
%         %      movefile([path, file{i}], ['C:\Users\Gama\Documents\Github\Morphology\all\', file{i}])
%     end
%     ini=j+ini;
%     
%     
% end


pathi='C:\Users\Gama\Documents\new_images\database\';
f=dir([pathi,'*.jpg']);
file={};
for i=1:length(f)
    file(i)={f(i).name};
end
ini=length(file);

pathi='C:\Users\Gama\Documents\new_images\database2\';
f=dir([pathi,'*.jpg']);
file={};
for i=1:length(f)
    file(i)={f(i).name};
end
path=pathi;

for i=1:length(indi)
    im = imread([path, file{indi(i)}]);

        imwrite(im, ['C:\Users\Gama\Documents\new_images\database\', num2str(i+ini),'.jpg']);
      %    movefile([path, file{i}], ['C:\Users\Gama\Documents\Github\Morphology\all\', file{i}])
end

    





