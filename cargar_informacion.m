armon = dir('armonicos*.mat'); 
[~, Index] = natsort({armon.name});
armon   = armon(Index);

excen = dir('excentricidad*.mat'); 
[~, Index] = natsort({excen.name});
excen   = excen(Index);

load('allclass/roundness.mat');

load('allclass/sphericity.mat'); 
roundness = roundness';
n = size(roundness,1);
N = 40; %Numero de armonicos

armonicosCell = zeros([n 4*N]);
cellExcentricity = zeros([n N]);
%roundClassifications = zeros([n,1]);
%spherClassifications = zeros([n,1]);



for q = 1:length(excen) 
    load(armon(q).name); 
    load(excen(q).name); 
    %load(round(q).name);
    %load(spher(q).name);
    armonicosCell((q-1)*50+1:q*50,:) = eval(sprintf('armonicos_database%d',q)); 
    cellExcentricity((q-1)*50+1:q*50,:) = eval(sprintf('excentricidad_database%d',q)); 
    %roundClassifications((q-1)*50+1:q*50) = roundness;
    %pherClassifications((q-1)*50+1:q*50) = eval(sprintf('sphericity_database%d',q));
end

finalCellRound = [num2cell(roundness) num2cell(armonicosCell)]; %Combinamos las celdas de excentricidad y clasificaciones
finalCellSpher = [num2cell(sphericity) num2cell(cellExcentricity)]; %Combinamos las celdas de excentricidad y clasificaciones
tableExcen = cell2table(finalCellSpher);
tableArmon = cell2table(finalCellRound);
tableExcen.Properties.VariableNames{1} = 'sphericity';
tableArmon.Properties.VariableNames{1} = 'roundness';
writetable(tableExcen,"sphericitydata.xlsx",'Sheet',1,'Range','A1')
writetable(tableArmon,"roundnessdata.xlsx",'Sheet',1,'Range','A1')