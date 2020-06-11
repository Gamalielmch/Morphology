function spherClassifications = exportar_esfericidad
    files = dir('allclass/*.jpg');  % specigy the extension of your image file
    [~, Index] = natsort({files.name});
    files   = files(Index);
    n = numel(files);
    spherClassifications = zeros([n,1]);
    for imnum = 351:700
        if (mod(imnum,25) == 0)
           imnum 
        end
        filename = strcat(files(imnum).folder,'\',files(imnum).name);
        spher = sphericity(filename,0);
        spherClassifications(imnum) = spher;
    end   
end