function spherClassifications = exportar_esfericidad
    files = dir('database14/*.jpg');  % specigy the extension of your image file
    [~, Index] = natsort({files.name});
    files   = files(Index);

    n = numel(files);
    spherClassifications = zeros([n,1]);
    
    for imnum = 1:numel(files)
        filename = strcat(files(imnum).folder,'\',files(imnum).name);
        spher = sphericity(filename,0);
        spherClassifications(imnum) = spher;
    end
        
end