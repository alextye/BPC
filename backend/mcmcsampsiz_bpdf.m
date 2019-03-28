% mcmcsampsiz - run mcmc such that effective sample size reaches target value
%
%
function [pchain,logLk,effSampN,accVec] = mcmcsampsiz_bpdf(NmcEff,logLkFun,...
        sig2,pchainTune,logLkTune,effSampNTune,accVecTune,method,areax,area_basis,fileID)

    sig2in = sig2;
    
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
    
    [pchainB,logLkB,accVecB]=mcmcfun(NchainB,pchainTune(end,:),logLkFun,sig2,areax,area_basis);
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

        [pchainRem,logLkRem,accVecRem]=...
            mcmcfun(NchainRem,pchain(end,:),logLkFun,sig2,areax,area_basis);
        
        %eliminate the first 10% of the chain.
        np = size(pchain,1);
        pchain = [pchain([ceil(np*.1):end],:); pchainRem];
        logLk = [logLk([ceil(np*.1):end],:); logLkRem];
        
        accVec = (accVec*(length(logLk)-NchainRem) + ...
                  accVecRem*NchainRem)/length(logLk);
        fprintf(fileID, '%s\n\n', strcat('accVec = ',mat2str(accVec)));
        [effSampN,dump,dump,dump,dump,nc_flag] = mcmcautocorr_bpdf(pchain);
        fprintf(fileID, '%s\n\n', strcat('effSampN = ',mat2str(effSampN)));
        fprintf(fileID, '%s\n\n', strcat('nc_flag = ',mat2str(nc_flag)));
%         while(nc_flag==1)
%             [pchainRem,logLkRem,accVecRem]=...
%                 mcmcfun(NchainRem,pchain(end,:),logLkFun,sig2,areax,area_basis);
%             pchain = [pchain;pchainRem];
%             logLk = [logLk; logLkRem];
%             accVec = (accVec*(length(logLk)-NchainRem) + ...
%                 accVecRem*NchainRem)/length(logLk)
%             [effSampN,dump,dump,dump,dump,nc_flag] = mcmcautocorr_bpdf(pchain);
%             effSampN
%             nc_flag
%         end

        
        NchainTarget = max(ceil(NmcEff*effSampN));
        fprintf(fileID, '%s\n\n', strcat('NchainTarget = ',mat2str(NchainTarget)));
        NchainRem = NchainTarget-size(pchain,1);
        fprintf(fileID, '%s\n\n', strcat('NchainRem = ',mat2str(NchainRem)));
        if nc_flag==1
            NchainRem = 0;
        end
    end
    fprintf(fileID, '%s\n\n', strcat('NmcEff = ',mat2str(NmcEff)));
    fprintf(fileID, '%s\n\n', strcat('effSampN = ',mat2str(effSampN)));
end