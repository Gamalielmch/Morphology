%function exportar_informacion
    pathi=[pwd,'\testimage'];
    [file,path] = uigetfile({'*.jpg;*.png;*.bmp;*.jpeg;'},'Select One or More Files','MultiSelect', 'on',pathi);
    if iscell(file)
        n=length(file);
    else
        n=1;
        file={file};
    end
    %images ='testimage';
    filename = 'excentricidad_data.xlsx'; 
    %jpgfiles = dir(fullfile(images,'\*.jpg'));
    
    %n = numel(jpgfiles); %Numero de imagenes en el directorio im_pruebas
    N = 40;

    cellExcentricity = zeros([n N*4]);
    roundClassifications = cell([n,1]);
    spherClassifications = cell([n,1]);
    %shper = zeros(1,n);
    
    %roundness=zeros(1,n);
    for imnum = 1:n
        spher = sphericity([path, file{imnum}],1);
        cons = armonicosFourierEliptico(N,[path, file{imnum}]);
        roundness = inscribedCircles_5([path, file{imnum}],N,15,file{imnum});
       
        rClass = roundness; %Se obtiene el valor del primer decimal, debido a que solo existen 9 clases 0.1 - 0.9
        roundClassifications{imnum} = rClass;
        spherClassifications{imnum} = spher;
        cellExcentricity(imnum,:) = cons';
    end
    
    finalCell = [roundClassifications spherClassifications num2cell(cellExcentricity)]; %Combinamos las celdas de excentricidad y clasificaciones
    tableExcen = cell2table(finalCell);
    writetable(tableExcen,filename,'Sheet',1,'Range','A1')
%end

