%SGD Method actualiza los pesos en cada entrada, por lo que los va a
%actualizar n veces, dependiendo de los n datos
function weight = SGD_method (weight, input, correct_output)
    alpha = 0.9; %learning rate
    N = size(input,1);
    for k = 1:N
        t_input = input(k,:)';
        d = correct_output(k);
        
        w_sum = weight * t_input;
        output = sigmoid(w_sum);
        error = d - output;
        delta = output * (1 - output) * error; %Delta rule
        
        
        dWeight = alpha * delta * t_input;
        
        weight = weight + dWeight';
    end
end