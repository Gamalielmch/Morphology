%function exportar_informacion
files = dir('database1/*.jpg');  % specigy the extension of your image file
[~, Index] = natsort({files.name});
files   = files(Index);

files2 = dir('database2/*.jpg');  % specigy the extension of your image file
[~, Index] = natsort({files2.name});
files2   = files2(Index);

files3 = dir('database3/*.jpg');  % specigy the extension of your image file
[~, Index] = natsort({files3.name});
files3   = files3(Index);

files_comb = [files;files2;files3];

r1 = load('roundness_database1.mat');
r2 = load('roundness_database2.mat');
r3 = load('roundness_database3.mat');
roundness = [r1.roundness,r2.roundness,r3.roundness];

% pathi=[pwd,'\database1'];
%     [file,path] = uigetfile({'*.jpg;*.png;*.bmp;*.jpeg;'},'Select One or More Files','MultiSelect', 'on',pathi);
%     if iscell(file)
%         n=length(file);
%     else
%         n=1;
%         file={file};
%     end
    %images ='testimage';
    filename = 'excentricidad_data.xlsx'; 
    %jpgfiles = dir(fullfile(images,'\*.jpg'));
    
    %n = numel(jpgfiles); %Numero de imagenes en el directorio im_pruebas
    N = 40;
    n = numel(files_comb);
    cellExcentricity = zeros([n N]);
    roundClassifications = cell([n,1]);
    spherClassifications = cell([n,1]);
    %shper = zeros(1,n);
    
    %roundness=zeros(1,n);
    for imnum = 1:numel(files_comb)
        filename = strcat(files_comb(imnum).folder,'\',files_comb(imnum).name);
       % image = imread(filename);
        spher = sphericity(filename,0);
        cons = armonicosFourierEliptico(N,filename);
        %roundness = inscribedCircles_5([path, file{imnum}],N,15,file{imnum});
       
        %rClass = roundness; %Se obtiene el valor del primer decimal, debido a que solo existen 9 clases 0.1 - 0.9
        %roundClassifications{imnum} = rClass;
        spherClassifications{imnum} = spher;
        cellExcentricity(imnum,:) = cons';
    end
    
    finalCell = [num2cell(roundness') spherClassifications num2cell(cellExcentricity)]; %Combinamos las celdas de excentricidad y clasificaciones
    tableExcen = cell2table(finalCell);
    writetable(tableExcen,filename,'Sheet',1,'Range','A1')
%end

