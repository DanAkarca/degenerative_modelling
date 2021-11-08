% the degenerative sptl model
% written by danyal akarca, university of cambridge, 2021
function b = fcn_sptl(A,D,mend,eta,modelvar)
% number of nodes
n = length(D);
% starting number of edges
mseed = nnz(A)/2;
% get the model variable
switch modelvar
    case 'powerlaw'
        Fd = D.^eta;
    case 'exponential'
        Fd = exp(eta*D);
end
% list all edges
[u,v] = find(triu(ones(n),1));
% get all indices
indx = (v - 1)*n + u;
% get the probability of all current edges
P = Fd(indx);
% initialise
b = zeros(mseed-mend,1);
% loop over the degenerative model
for i = 1:(mseed-mend);
    % form the cumulative probability 
    C = [0; cumsum(P)];
    % pick an edge according to this probability
    r = sum(rand*C(end) >= C);
    % take it away from the network
    A(u(r),v(r)) = 0;
    A(v(r),u(r)) = 0;
    % update probability
    P(r) = 0;
    % keep in order as they were removed so we can recapitulate developmental trends
    C = zeros(n); C(u(r),v(r))=1;
    b(i) = find(C);
    % display
    disp(sprintf('sptl model edge %g of %g removed',i,mseed-mend));
end
end