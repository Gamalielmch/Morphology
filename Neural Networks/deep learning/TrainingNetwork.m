%Deep neural network de 3 capas ocultas y 1 de salida para reconocer
%numeros en matrices binarias
input_Image = zeros(5,5,5);

input_Image(:,:,1) = [1 0 0 1 1;
                      1 1 0 1 1;
                      1 1 0 1 1;
                      1 1 0 1 1;
                      1 0 0 0 1];
           
input_Image(:,:,2) = [0 0 0 0 1;
                      1 1 1 1 0;
                      1 0 0 0 1;
                      0 1 1 1 1;
                      0 0 0 0 0];
                
input_Image(:,:,3) = [0 0 0 0 1;
                      1 1 1 1 0;
                      1 0 0 0 1;
                      1 1 1 1 0;
                      0 0 0 0 0];
                  
input_Image(:,:,4) = [1 1 1 0 1;
                      1 1 0 0 1;
                      1 0 1 0 1;
                      0 0 0 0 0;
                      1 1 1 0 1];
                  
input_Image(:,:,5) = [0 0 0 0 0;
                      0 1 1 1 1;
                      0 0 0 0 1;
                      1 1 1 1 0;
                      0 0 0 0 1];
                  
correct_input = [1 0 0 0 0;
                 0 1 0 0 0;
                 0 0 1 0 0;
                 0 0 0 1 0;
                 0 0 0 0 1];
 %primera capa tiene 25 nodos,
 %segunda capa tiene 20
 %tercera capa tiene 20
 %capa de salida tiene 5 nodos, porque son las respectivas salidas que puede tomar
w = {2*rand(20,25)-1, 2*rand(20,20)-1, 2*rand(20,20)-1,2*rand(5,20)-1};

for epoch = 1:10000
   w = DeepLearning(w,input_Image,correct_input); 
end

save('DeepNeuralNetwork.mat');