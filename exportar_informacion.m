function exportar_informacion
    
    images ='testimage';
    filename = 'excentricidad_data.xlsx'; 
    jpgfiles = dir(fullfile(images,'\*.jpg'));
    
    n = numel(jpgfiles); %Numero de imagenes en el directorio im_pruebas
    N = 45;

    cellExcentricity = zeros([n N]);
    roundClassifications = cell([n,1]);
    
    for imnum = 1:n
        imname = jpgfiles(imnum).name;
        img_route = fullfile(images,imname);
       % excen = armonicosFourierEliptico(N,img_route);
        roundness = inscribedCircles_3(img_route,N,.1,20,imname);
        
        rClass = fix(roundness*10)/10; %Se obtiene el valor del primer decimal, debido a que solo existen 9 clases 0.1 - 0.9
        roundClassifications{imnum} = rClass;
     %   cellExcentricity(imnum,:) = excen';
    end
    
    finalCell = [classifications  num2cell(cellExcentricity)]; %Combinamos las celdas de excentricidad y clasificaciones
    tableExcen = cell2table(finalCell);
    writetable(tableExcen,filename,'Sheet',1,'Range','A1')
end

