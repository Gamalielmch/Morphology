function [armonClassifications,excenClassifications] = exportar_armonicos
    files = dir('allclass/*.jpg');  % specigy the extension of your image file
    [~, Index] = natsort({files.name});
    files   = files(Index);

    n = numel(files);
    N = 40;
    armonClassifications = zeros([n,N*4]);
    excenClassifications = zeros([n,N]);
    
    for imnum = 1:numel(files)
        if (mod(imnum,25) == 0)
           imnum 
        end
        filename = strcat(files(imnum).folder,'\',files(imnum).name);
        [cons,excen] = armonicosFourierEliptico(N,filename);
        armonClassifications(imnum,:) = cons;
        excenClassifications(imnum,:) = excen;
    end
        
end