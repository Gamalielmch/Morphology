load('Trained_network.mat');
input = [0 0 1;
         0 1 1;
         1 0 1;
         1 1 1];
N = size(input,1);

for k = 1:N
   t_input = input(k,:)';
   weighted_sum = weight * t_input;
   ouput = sigmoid(weighted_sum)
end