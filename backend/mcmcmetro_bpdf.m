%mcmcmetro- runs a mcmc with metropolis-hastings using full covariance matrix
%
%create mcmc chain of length Nchain using metropolis-hastings algorithm with
%multivariate normal proposal distribution.
%
%NOTE: the chain stores resulting iterations (so it does not store the initial
%value, unless the value is repeated). This is so that the function can be
%called repeatedly using the last iteration, and all the chains concatenated
%simply to form one long valid chain.
function [pchain,logLk,accRatio] = mcmcmetro_bpdf(Nchain,pinit,logLkFun,sig2,areax,area_basis)

    pinit = pinit(:)';
    if(size(sig2,1)~=size(sig2,2))
        sig2 = diag(sig2); 
    end

    Np = length(pinit);
    assert(size(sig2,1)==Np,...
        'length of parameter array must equal length of sig2 array');

    pchain = nan*ones(Nchain,Np);
    logLk  = nan*ones(Nchain,1);
    accHist= nan*ones(Nchain,1);

    EPSarea = 1e-5;
    larea = sp_log_area(pinit,areax,area_basis);
    while(abs(larea)>EPSarea)
        pinit = pinit - larea;
        larea = sp_log_area(pinit,areax,area_basis);
    end
    
    prev = pinit;
%    csvwrite('temp_dump/pinit',pinit);
    prevLogLk = -1*logLkFun(pinit);
    logmaxp = prevLogLk;
  %  propparam = ones(Nchain,Np);
  %  propchain = ones(Nchain,1);

    for(i=1:Nchain)
        prop = mvnrnd(prev,sig2,1);
        prop(find(prop<-50)) = -50;
        prop = prop - sp_log_area(prop,areax,area_basis);
        
        %propparam(i,:) = prop;
        propLogLk = -1*logLkFun(prop);
        if(~isnan(propLogLk))
        %    break;
            paccept = min(1,exp(propLogLk-prevLogLk));
        else
            paccept = 0;
        end
        %propchain(i) = propLogLk;
        
        %could use the maximum log_prob_dens value instead of 1?
        
        if(rand(1,1) < paccept)
            % store prop value in chain
            pchain(i,:) = prop;
            logLk(i)    = propLogLk;
            accHist(i)  = true;

            % update prev values
            prev      = prop;
            prevLogLk = propLogLk;
        else
            % repeat prev value in chain
            pchain(i,:) = prev;
            logLk(i)    = prevLogLk;
            accHist(i)  = false;

            % do not change prev values
        end
        %if(mod(i, 100)==0)
            %i
        %end
    end
    accRatio = mean(accHist);
  %  csvwrite('temp_dump/pchain',pchain);
  %  csvwrite('temp_dump/propchain',propchain);
  %  csvwrite('temp_dump/propparam',propparam);
  %  csvwrite('temp_dump/mcmcsig2',sig2);
end