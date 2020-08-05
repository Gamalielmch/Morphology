net=feedforwardnet([20 20 20 20 20 20 20]);
%armonics = armonics{:,:};
%roundnessdata = roundnessdata{:,:};
armon = armonics(:,41:80);
armon = zscore(armon);
coeff = pca(armon);
armon = armon * coeff;
[net,tr] = train(net,armon',roundnessdata');
net.layers{1}.transferFcn = 'poslin';
net.layers{2}.transferFcn = 'poslin';
net.layers{3}.transferFcn = 'poslin';
net.layers{4}.transferFcn = 'poslin';
net.layers{5}.transferFcn = 'poslin';
net.layers{6}.transferFcn = 'poslin';
net.layers{7}.transferFcn = 'poslin';
outputs = net(armon')
%view(net);
% Plots
% Uncomment these lines to enable various plots.
% figure, plotsomtop(net)
% figure, plotsomnc(net)
% figure, plotsomnd(net)
% figure, plotsomplanes(net)
% figure, plotsomhits(net,inputs)
% figure, plotsompos(net,inputs)