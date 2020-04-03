input = [0 0 1;
         0 1 1;
         1 0 1;
         1 1 1];
correct_output = [0 0 1 1];
weight = 2*rand(1,3) - 1;

for epoch = 1:10000
   weight = SGD_method(weight, input, correct_output); 
end

save('Trained_network.mat');