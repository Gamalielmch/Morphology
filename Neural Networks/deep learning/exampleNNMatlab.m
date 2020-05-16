net = feedforwardnet;
net.numInputs = 25;
net.numLayers = 4;
input_Image = zeros(5,25);

input_Image(1,:) = [1 0 0 1 1 1 1 0 1 1 1 1 0 1 1 1 1 0 1 1 1 0 0 0 1];
           
input_Image(2,:) = [0 0 0 0 1 1 1 1 1 0 1 0 0 0 1 0 1 1 1 1 0 0 0 0 0];
                
input_Image(3,:) = [0 0 0 0 1 1 1 1 1 0 1 0 0 0 1 1 1 1 1 0 0 0 0 0 0];
                  
input_Image(4,:) = [1 1 1 0 1 1 1 0 0 1  1 0 1 0 1 0 0 0 0 0 1 1 1 0 1];
                  
input_Image(5,:) = [0 0 0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1 1 0 0 0 0 0 1];
                  
correct_input = [1 1 1 1 1];

