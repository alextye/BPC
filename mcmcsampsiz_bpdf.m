% mcmcsampsiz - run mcmc such that effective sample size reaches target value
%
%
function [pchain,logLk,effSampN,accVec] = mcmcsampsiz_bpdf(NmcEff,logLkFun,...
        sig2,pchainTune,logLkTune,effSampNTune,accVecTune,method,knots,N_pts,fileID)

    sig2in = sig2;
    %factor by which to modify the step size--we want the accuracy vector
    %to approach 0.24.
%    sigfactor = accVecTune/0.24;
    
    if(strcmp(method,'gibbs'))
        mcmcfun = @mcmcgibbs;
    elseif(strcmp(method,'metro'))
        mcmcfun = @mcmcmetro_bpdf;
    else
        disp('Error in mcmcsampsiz:');
        assert(false,'method must be ''gibbs'' or ''metro''.');
    end

    NchainTarget = max(ceil(NmcEff*effSampNTune));
    fprintf(fileID, '%s\n\n', strcat('NchainTarget = ',mat2str(NchainTarget)));
    NchainTune = size(pchainTune,1);

    if(~isreal(NchainTarget))
    keyboard;
    end
    
    NchainB = ceil(0.5*(NchainTarget-NchainTune));
    fprintf(fileID, '%s\n\n', strcat('NchainB = ',mat2str(NchainB)));
%    keyboard;
    
    [pchainB,logLkB,accVecB]=mcmcfun(NchainB,pchainTune(end,:),logLkFun,sig2,knots,N_pts);
    pchain = [pchainTune; pchainB];
    logLk = [logLkTune; logLkB];
    accVec = (accVecTune*length(logLkTune) + ...
        accVecB*length(logLkB))/length(logLk);
    fprintf(fileID, '%s\n\n', strcat('accVec = ',mat2str(accVec)));
    [effSampN] = mcmcautocorr_bpdf(pchain);
    fprintf(fileID, '%s\n\n', strcat('effSampN = ',mat2str(effSampN)));
    NchainTarget = max(ceil(NmcEff*effSampN));
    
    if(~isreal(NchainTarget))
    keyboard;
    end

    NchainRem = NchainTarget-(NchainTune+NchainB);
    fprintf(fileID, '%s\n\n', strcat('NchainRem = ',mat2str(NchainRem)));

    while(NchainRem > 0)
%        keyboard;
%        sig2 = sig2in;
%        sigfactor = accVec/0.24;
        [pchainRem,logLkRem,accVecRem]=...
            mcmcfun(NchainRem,pchain(end,:),logLkFun,sig2,knots,N_pts);
        pchain = [pchain;pchainRem];
        logLk = [logLk; logLkRem];
        accVec = (accVec*(length(logLk)-NchainRem) + ...
                  accVecRem*NchainRem)/length(logLk);
        fprintf(fileID, '%s\n\n', strcat('accVec = ',mat2str(accVec)));
        [effSampN,dump,dump,dump,dump,nc_flag] = mcmcautocorr_bpdf(pchain);
        fprintf(fileID, '%s\n\n', strcat('effSampN = ',mat2str(effSampN)));
        fprintf(fileID, '%s\n\n', strcat('nc_flag = ',mat2str(nc_flag)));
        while(nc_flag==1)
            [pchainRem,logLkRem,accVecRem]=...
                mcmcfun(NchainRem,pchain(end,:),logLkFun,sig2,knots,N_pts);
            pchain = [pchain;pchainRem];
            logLk = [logLk; logLkRem];
            accVec = (accVec*(length(logLk)-NchainRem) + ...
                accVecRem*NchainRem)/length(logLk)
            [effSampN,dump,dump,dump,dump,nc_flag] = mcmcautocorr_bpdf(pchain);
            effSampN
            nc_flag
        end
        NchainTarget = max(ceil(NmcEff*effSampN));
        fprintf(fileID, '%s\n\n', strcat('NchainTarget = ',mat2str(NchainTarget)));
        NchainRem = NchainTarget-size(pchain,1);
        fprintf(fileID, '%s\n\n', strcat('NchainRem = ',mat2str(NchainRem)));
        %keyboard;
    end
    fprintf(fileID, '%s\n\n', strcat('NmcEff = ',mat2str(NmcEff)));
    fprintf(fileID, '%s\n\n', strcat('effSampN = ',mat2str(effSampN)));
end