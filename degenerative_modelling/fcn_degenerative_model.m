% the degenerative sptl and value degenerative model
% written by danyal akarca, university of cambridge, 2021
function b = fcn_degenerative_model(A,K,D,model,mend,eta,gam,modelvar,epsilon)
% take values of the seed
K = K + epsilon;
% number of nodes
n = length(D);
% number of starting edges
mseed = nnz(A)/2;
% take the model variables
mv1 = modelvar{1};
mv2 = modelvar{2};
switch mv1
    case 'powerlaw'
        Fd = D.^eta;
    case 'exponential'
        Fd = exp(eta*D);
end
switch mv2
    case 'powerlaw'
        Fk = K.^gam;
    case 'exponential'
        Fk = exp(gam*K);
end
% probabilty functions
Ff = Fd.*Fk;
[u,v] = find(triu(ones(n),1));
indx = (v - 1)*n + u;
% probability
P = Ff(indx);
% initialise
b = zeros(mseed-mend,1);
% loop over edges to remove
for i = 1:(mseed-mend);
    % sample according to the probability distribution
    e = randsample(indx,1,true,P);
    r = find(indx==e);
    % remove edges
    uu = u(r);
    vv = v(r);
    A(uu,vv) = 0;
    A(vv,uu) = 0;
    % update the matrix depending on the model
    switch model
        case 'matching'
            K = matching_ind(A);
            K = K + K';
        case 'neighbours'
            K = (A*A)*~eye(n);
        case 'clu-avg'
            K = clustering_coef_bu(A);
            [ii jj] = meshgrid(K,K);
            S = mean([ii(:),jj(:)],2);
            K = reshape(S,[n n]);
        case 'clu-min'
            K = clustering_coef_bu(A);
            [ii jj] = meshgrid(K,K);
            S = min([ii(:),jj(:)],2);
            K = reshape(S,[n n]);
        case 'clu-max'
            K = clustering_coef_bu(A);
            [ii jj] = meshgrid(K,K);
            S = max([ii(:),jj(:)],2);
            K = reshape(S,[n n]);
        case 'clu-diff'
            K = clustering_coef_bu(A);
            [ii jj] = meshgrid(K,K);
            S = abs([ii(:)-jj(:)]);
            K = reshape(S,[n n]);
        case 'clu-prod'
            K = clustering_coef_bu(A);
            [ii jj] = meshgrid(K,K);
            S = [ii(:).*jj(:)]
            K = reshape(S,[n n]);
        case 'deg-avg'
            K = degrees_und(A);
            [ii jj] = meshgrid(K,K);
            S = mean([ii(:),jj(:)],2);
            K = reshape(S,[n n]);
        case 'deg-min'
            K = degrees_und(A);
            [ii jj] = meshgrid(K,K);
            S = min([ii(:),jj(:)],2);
            K = reshape(S,[n n]);
        case 'deg-max'
            K = degrees_und(A);
            [ii jj] = meshgrid(K,K);
            S = max([ii(:),jj(:)],2);
            K = reshape(S,[n n]);
        case 'deg-diff'
            K = degrees_und(A);
            [ii jj] = meshgrid(K,K);
            S = abs([ii(:)-jj(:)]);
            K = reshape(S,[n n]);
        case 'deg-prod'
            K = degrees_und(A);
            [ii jj] = meshgrid(K,K);
            S = [ii(:).*jj(:)]
            K = reshape(S,[n n]);
    end
    switch mv2
        case 'powerlaw'
            Fk = K.^gam;
        case 'exponential'
            Fk = exp(gam*K);
    end
    % update the probability
    P = Fd.*Fk.*A;
    P = P(indx);
    % remove nans if they arise (does happen sometimes)
    P(isnan(P))=0;
    % keep in order as they were removed so we can recapitulate developmental trends
    C = zeros(n); C(uu,vv)=1;
    b(i) = find(C);
    % display
    disp(sprintf('%s degenerative model %g edges remain',model,nnz(A)/2));
end
end