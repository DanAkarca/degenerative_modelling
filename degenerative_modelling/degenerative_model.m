function b = degenerative_model(A,D,mend,modeltype,modelvar,params,epsilon)
%{
DEGENERATIVE_MODEL          Run generative model code
B = DEGENERATIVE_MODEL(A,D,mend,modeltype,modelvar,params,epsilon)

Generates synthetic networks using the models described in the study by
Akarca et al (2022) in (perhaps a future memory of being in a published journal).
Inputs:
           A,          starting binary network of connections
           D,          cost matrix (e.g. Euclidean distance/fiber length matrix)
           mend,          number of connections that should be present in
                       final synthetic network
           modeltype,  specifies the degenerative rule (see below)
           modelvar,   specifies whether the degenerative rules are based on
                       power-law or exponential relationship
                       ({'powerlaw'}|{'exponential})
           params,     either a vector (in the case of the geometric
                       model) or a matrix (for all other models) of
                       parameters at which the model should be evaluated.
           epsilon,    the baseline probability of removing a particular
                       connection (should be a very small number
                       {default = 1e-5}).

Output:
           B,          number of networks matrix of connections x m


Full list of model types:
(each model type realizes a different degenerative rule)

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

Reference: Akarca et al (2022) (future memory of being in a published journal)
Adapted from generative_model.m written by Richard Betzel, 2015
%}
% if epislon hasn't been set, set it
if ~exist('epsilon','var')
    epsilon = 1e-5;
end
% compute the number of nodes
n = length(D);
% compute the number of starting edges
mstart = nnz(A)/2;
% compute the number of parameters
nparams = size(params,1);
% initialise then output array
b = zeros(nparams,mstart-mend);
% swtich modeltype and run the degenerative model
switch modeltype
    case 'matching'
        % compute the initial seed values
        Kseed = matching_ind(A);
        Kseed = Kseed + Kseed';
        % loop over parameters
        for iparam = 1:nparams
            % take eta and gamma
            eta = params(iparam,1);
            gam = params(iparam,2);
            % run the matching degenerative model
            b(iparam,:) = fcn_degenerative_model(A,Kseed,D,'matching',mend,eta,gam,modelvar,epsilon);
            % display
            disp(sprintf('computed %s parameter combination %g of %g',modeltype,iparam,nparams));
        end
    case 'neighbours'
        % compute the initial seed values
        Kseed = (A*A).*~eye(n);
        % loop over parameters
        for iparam = 1:nparams
            % take eta and gamma
            eta = params(iparam,1);
            gam = params(iparam,2);
            % run the matching degenerative model
            b(iparam,:) = fcn_degenerative_model(A,Kseed,D,'neighbours',mend,eta,gam,modelvar,epsilon);
            % display
            disp(sprintf('computed %s parameter combination %g of %g',modeltype,iparam,nparams));
        end
    case 'clu-avg'
        % compute the initial seed values
        Kseed = clustering_coef_bu(A);
        [ii jj] = meshgrid(Kseed,Kseed);
        S = mean([ii(:),jj(:)],2);
        Kseed = reshape(S,[n n]);
        % loop over parameters
        for iparam = 1:nparams
            % take eta and gamma
            eta = params(iparam,1);
            gam = params(iparam,2);
            % run the matching degenerative model
            b(iparam,:) = fcn_degenerative_model(A,Kseed,D,'clu-avg',mend,eta,gam,modelvar,epsilon);
            % display
            disp(sprintf('computed %s parameter combination %g of %g',modeltype,iparam,nparams));
        end
    case 'clu-min'
        % compute the initial seed values
        Kseed = clustering_coef_bu(A);
        [ii jj] = meshgrid(Kseed,Kseed);
        S = min([ii(:),jj(:)],2);
        Kseed = reshape(S,[n n]);
        % loop over parameters
        for iparam = 1:nparams
            % take eta and gamma
            eta = params(iparam,1);
            gam = params(iparam,2);
            % run the matching degenerative model
            b(iparam,:) = fcn_degenerative_model(A,Kseed,D,'clu-min',mend,eta,gam,modelvar,epsilon);
            % display
            disp(sprintf('computed %s parameter combination %g of %g',modeltype,iparam,nparams));
        end
    case 'clu-max'
        % compute the initial seed values
        Kseed = clustering_coef_bu(A);
        [ii jj] = meshgrid(Kseed,Kseed);
        S = max([ii(:),jj(:)],2);
        Kseed = reshape(S,[n n]);
        % loop over parameters
        for iparam = 1:nparams
            % take eta and gamma
            eta = params(iparam,1);
            gam = params(iparam,2);
            % run the matching degenerative model
            b(iparam,:) = fcn_degenerative_model(A,Kseed,D,'clu-max',mend,eta,gam,modelvar,epsilon);
            % display
            disp(sprintf('computed %s parameter combination %g of %g',modeltype,iparam,nparams));
        end
    case 'clu-diff'
        % compute the initial seed values
        Kseed = clustering_coef_bu(A);
        [ii jj] = meshgrid(Kseed,Kseed);
        S = abs(ii(:)-jj(:));
        Kseed = reshape(S,[n n]);
        % loop over parameters
        for iparam = 1:nparams
            % take eta and gamma
            eta = params(iparam,1);
            gam = params(iparam,2);
            % run the matching degenerative model
            b(iparam,:) = fcn_degenerative_model(A,Kseed,D,'clu-diff',mend,eta,gam,modelvar,epsilon);
            % display
            disp(sprintf('computed %s parameter combination %g of %g',modeltype,iparam,nparams));
        end
    case 'clu-prod'
        % compute the initial seed values
        Kseed = clustering_coef_bu(A);
        [ii jj] = meshgrid(Kseed,Kseed);
        S = ii(:).*jj(:);
        Kseed = reshape(S,[n n]);
        % loop over parameters
        for iparam = 1:nparams
            % take eta and gamma
            eta = params(iparam,1);
            gam = params(iparam,2);
            % run the matching degenerative model
            b(iparam,:) = fcn_degenerative_model(A,Kseed,D,'clu-diff',mend,eta,gam,modelvar,epsilon);
            % display
            disp(sprintf('computed %s parameter combination %g of %g',modeltype,iparam,nparams));
        end
    case 'deg-avg'
        % compute the initial seed values
        Kseed = degrees_und(A);
        [ii jj] = meshgrid(Kseed,Kseed);
        S = mean([ii(:),jj(:)],2);
        Kseed = reshape(S,[n n]);
        % loop over parameters
        for iparam = 1:nparams
            % take eta and gamma
            eta = params(iparam,1);
            gam = params(iparam,2);
            % run the matching degenerative model
            b(iparam,:) = fcn_degenerative_model(A,Kseed,D,'deg-avg',mend,eta,gam,modelvar,epsilon);
            % display
            disp(sprintf('computed %s parameter combination %g of %g',modeltype,iparam,nparams));
        end
    case 'deg-min'
        % compute the initial seed values
        Kseed = degrees_und(A);
        [ii jj] = meshgrid(Kseed,Kseed);
        S = min([ii(:),jj(:)],2);
        Kseed = reshape(S,[n n]);
        % loop over parameters
        for iparam = 1:nparams
            % take eta and gamma
            eta = params(iparam,1);
            gam = params(iparam,2);
            % run the matching degenerative model
            b(iparam,:) = fcn_degenerative_model(A,Kseed,D,'deg-min',mend,eta,gam,modelvar,epsilon);
            % display
            disp(sprintf('computed %s parameter combination %g of %g',modeltype,iparam,nparams));
        end
    case 'deg-max'
        % compute the initial seed values
        Kseed = degrees_und(A);
        [ii jj] = meshgrid(Kseed,Kseed);
        S = max([ii(:),jj(:)],2);
        Kseed = reshape(S,[n n]);
        % loop over parameters
        for iparam = 1:nparams
            % take eta and gamma
            eta = params(iparam,1);
            gam = params(iparam,2);
            % run the matching degenerative model
            b(iparam,:) = fcn_degenerative_model(A,Kseed,D,'deg-max',mend,eta,gam,modelvar,epsilon);
            % display
            disp(sprintf('computed %s parameter combination %g of %g',modeltype,iparam,nparams));
        end
    case 'deg-diff'
        % compute the initial seed values
        Kseed = degrees_und(A);
        [ii jj] = meshgrid(Kseed,Kseed);
        S = abs(ii(:)-jj(:));
        Kseed = reshape(S,[n n]);
        % loop over parameters
        for iparam = 1:nparams
            % take eta and gamma
            eta = params(iparam,1);
            gam = params(iparam,2);
            % run the matching degenerative model
            b(iparam,:) = fcn_degenerative_model(A,Kseed,D,'deg-diff',mend,eta,gam,modelvar,epsilon);
            % display
            disp(sprintf('computed %s parameter combination %g of %g',modeltype,iparam,nparams));
        end
    case 'deg-prod'
        % compute the initial seed values
        Kseed = degrees_und(A);
        [ii jj] = meshgrid(Kseed,Kseed);
        S = ii(:).*jj(:);
        Kseed = reshape(S,[n n]);
        % loop over parameters
        for iparam = 1:nparams
            % take eta and gamma
            eta = params(iparam,1);
            gam = params(iparam,2);
            % run the matching degenerative model
            b(iparam,:) = fcn_degenerative_model(A,Kseed,D,'deg-diff',mend,eta,gam,modelvar,epsilon);
            % display
            disp(sprintf('computed %s parameter combination %g of %g',modeltype,iparam,nparams));
        end    
    case 'sptl'
        for iparam = 1:nparams
            % take eta
            eta = params(iparam,1);
            % run the sptl degenerative model
            b(iparam,:) = fcn_degenerative_sptl_model(A,D,mend,eta,modelvar{1});
            % display
            disp(sprintf('computed %s eta %g of %g',modeltype,iparam,nparams));
        end
end
end