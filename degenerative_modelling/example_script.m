%% example script to run the degenerative_model
% written by danyal akarca, university of cambridge, 2021
% outputs the removed edges, in order of how they were removed
clear; clc;
%% addpath of your brain connectivity toolbox
bct_path = '/Users/da04/Desktop/PhD/Research/Toolboxes/BCT/2019_03_03_BCT/';
addpath(bct_path);
%% addpath of your degenerative modelling toolbox
degenerative_path = '/Users/da04/Desktop/degenerative_modelling/';
addpath(degenerative_path);
%% input an example starting network and cost term
% set the number of nodes (square)
nnode = 36;
% set the starting network
Astart = ones(nnode); Astart(find(eye(size(Astart,1))))=0;
% set the coordinates and cost function
C = [1:sqrt(nnode)]; [a b] = meshgrid(C,C); C = [a(:) b(:)];
D = squareform(pdist(C));
% compute initial number of connections
mstart = nnz(Astart)/2;
% set the parameters
params = [3 -0.5];
%% visualise the starting network
h = figure; h.Position = [100 100 900 250];
sgtitle('generative model hyperparameters');
subplot(1,3,1); 
imagesc(Astart);
b = gca; b.TickDir = 'out';
title('starting network adjacency');
subplot(1,3,2); 
imagesc(D); 
b = gca; b.TickDir = 'out';
title('cost');
subplot(1,3,3);
plot(graph(Astart),'XData',C(:,1),'YData',C(:,2),'NodeLabel',[]); 
b = gca; b.TickDir = 'out';
title('embedded network');
%% set the generative model hyperparameters
%{
set the model type
1.  'sptl'          spatial model
2.  'neighbors'     number of common neighbors
3.  'matching'      matching index
4.  'clu-avg'       average clustering coeff.
5.  'clu-min'       minimum clustering coeff.
6.  'clu-max'       maximum clustering coeff.
7.  'clu-diff'      difference in clustering coeff.
8.  'clu-prod'      product of clustering coeff.
9.  'deg-avg'       average degree
10. 'deg-min'       minimum degree
11. 'deg-max'       maximum degree
12. 'deg-diff'      difference in degree
13. 'deg-prod'      product of degree
%}
modeltype = 'matching';
% ending number of connections
mend = 100;
% set the model variables
modelvar = {'powerlaw','powerlaw'};
% set epsilon
epsilon = 1e-5;
%% run the degenerative model
b = degenerative_model(Astart,D,mend,modeltype,modelvar,params);
%% visualise simulation outputs
% select parameter to view
iparam = 1;
% form an empty network
Asynth = eye(nnode);
% populate removed edges in upper triangle
Asynth(b(iparam,:)) = 1;
% form complete network
Asynth = Asynth + Asynth';
% reverse to show the remaining network
Asynth = double(Asynth==0);
% plot
h = figure; h.Position = [100 100 1200 250];
subplot(1,4,1);
imagesc(Astart);
title('starting network');
subplot(1,4,2); 
plot(graph(Astart),'XData',C(:,1),'YData',C(:,2),'NodeLabel',[]);
subplot(1,4,3);
imagesc(Asynth);
sgtitle(sprintf('%s degenerative model | eta=%g, gamma=%g',modeltype,params(iparam,1),params(iparam,2)));
subplot(1,4,4); 
plot(graph(Asynth),'XData',C(:,1),'YData',C(:,2),'NodeLabel',[]);
%% visualise the trajectory of developmental change
% select parameter to view
iparam = 1;
% initialise statistics
nmeasures = 3;
statistics = zeros(mstart-mend,nmeasures);
statistics_labels = string({'modularity',...
    'efficiency','efficiency per weight'});
% loop over time
for t = 1:(mstart-mend)
    % form an empty network
    Asynth = eye(nnode);
    % populate removed edges in upper triangle
    Asynth(b(iparam,1:t)) = 1;
    % form complete network
    Asynth = Asynth + Asynth';
    % reverse to show the remaining network
    Asynth = double(Asynth==0);
    % plot statistics
    [~,statistics(t,1)] = modularity_und(Asynth);
    statistics(t,2) = efficiency_bin(Asynth,0);
    statistics(t,3) = statistics(t,2)/(nnz(Asynth)/2);
end
% initialise
h = figure; h.Position = [100 100 1000 250];
% plot
for stat = 1:nmeasures;
    subplot(1,nmeasures,stat);
    plot(statistics(:,stat),'linewidth',6); 
    ylabel(statistics_labels(stat));
    xlabel('time'); 
    k = gca; k.TickDir = 'out'; k.FontSize = 12;
    grid on;
end
sgtitle('global statistics during development');