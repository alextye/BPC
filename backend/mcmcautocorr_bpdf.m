% mcmcautocorr - estimate effective sample num based on mcmc autocorrelation
% 
% return Fisher transformation of correlation coefficient as logEffSampN
function [effSampN,effSampNErr,logEffSampN,logEffSampNErr,AR1coeff, nc_flag, fedup] = ...
        mcmcautocorr(pchain)

    fedup = 0;
    nc_flag = 0;
    Nchain = size(pchain,1);
    Nparam = size(pchain,2);

    %p = [z,mu,sig]
    %p(1) = z
    %p(2) = mu
    %p(3) = sig
    %Xt = pchain
%     negLogLkFun = @(p,Xt)(0.5*sum((Xt(2:end)-p(2)-...
%         (exp(2*p(1))-1)./(exp(2*p(1))+1)*(Xt(1:end-1)-p(2))).^2/exp(p(3))^2+...
%         2*p(3)));
% 
     mu = mean(pchain,1);
     sig = std(pchain-repmat(mu,Nchain,1));

%     AR1coeff = zeros(1,Nparam);
%     effSampN = zeros(1,Nparam);
%     effSampNErr = zeros(1,Nparam);
%     logEffSampN = zeros(1,Nparam);
%     logEffSampNErr = zeros(1,Nparam);
%     tic
%     for(i=1:Nparam)
%         dXt1 = pchain(2:end,i)-mu(i);
%         dXt0 = pchain(1:end-1,i)-mu(i);
%         iAR1 = sum(dXt1.*dXt0)./sum(dXt0.^2);
%         
%         if(iAR1>0.9999)
%             iAR1 = 0.95;
%             nc_flag = 1;
%         end
%         
%         z0 = 0.5*log((1+iAR1)./(1-iAR1));
%         p0 = [z0,mu(i),log(sig(i)*sqrt(1-iAR1^2))];
%         try
%             inegLogLkFun=@(p)(negLogLkFun(p,pchain(:,i))-negLogLkFun(p0,pchain(:,i)));
%             options = optimset('TolX',1e-8,'TolFun',1e-8,'LargeScale','off',...
%                 'Display','off');
%             [pfit fval flg output grd negHess] = fminunc(inegLogLkFun,p0,options);
%             %[negHess] = hessian(inegLogLkFun,pfit);
%             sig2 = inv(negHess);
%         catch
%             pfit = p0;
%             sig2 = NaN;
%         end 
%         AR1coeff(i) = (exp(2*pfit(1))-1)./(exp(2*pfit(1))+1);
%         %(1-AR1coeff(i))/(1+AR1coeff(i))*Nchain
%         effSampN(i) = exp(2*pfit(1));
%         
%         if(iAR1==0.95)
%             effSampN(i) = 1;
%         end
%         
%         if(~isreal(exp(2*pfit(1)))|exp(2*pfit(1))>1000000)
%             fedup = 1;
%         end
%         effSampNErr(i) = sqrt(sig2(1,1)*(2*exp(2*pfit(1)))^2);
%         logEffSampN(i) = pfit(1);
%         logEffSampNErr(i) = sqrt(sig2(1,1));
%         
%     end
%     toc
%    keyboard;
    
%    find MLE (max likelihood estimate) of autocorrelation coefficient for
%    each variable
    AR1coeff = zeros(1,Nparam);
    AR1err = zeros(1,Nparam);
    
    %tic
    
    for(i=1:Nparam)
       dXt1 = pchain(2:end,i)-mu(i);
       dXt0 = pchain(1:end-1,i)-mu(i);
       AR1coeff(i) = sum(dXt1.*dXt0)./sum(dXt0.^2);
       iHess = -sum(dXt0.^2)./sig(i)^2/(1-AR1coeff(i)^2);
       AR1err(i) = sqrt(-inv(iHess));
    end
    

    %effSampFac = (1-AR1coeff)./(1+AR1coeff);
    effSampN = (1+AR1coeff)./(1-AR1coeff);
    effSampN(find(AR1coeff>0.9999)) = 1;
    if length(find(AR1coeff>0.9999))>0
        nc_flag = 1;
    end
    effSampNErr = 2*AR1err./(1+AR1coeff).^2;
    logEffSampN = log(effSampN);
    logEffSampNErr = log(effSampNErr);
    
    %toc
end