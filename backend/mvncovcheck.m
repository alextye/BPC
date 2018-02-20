function [y, Dcheck] = mvncovcheck(opt, covmat, logLkfn)
%refines an estimated covariance matrix from the hessian from fminunc
%or similar function.  Iterates through eigenvectors and uses bisection
%algorithm to find more accurate eigenvalues than those input.

%PARAMETERS
%opt is the maximum likelihood model for a given zircon sample.
%covmat is the covariance matrix derived from the Hessian returned by the
%   optimization algorithm.
%logLkfn is the handle of the function used to calculate the likelihoods of
%   different models.

%OUTPUT
%y is the covariance matrix, in the original basis, after the correction
%procedure.
%Dcheck is the set of updated eigenvalues.

try
    %find eigenvectors and eigenvalues of the covariance matrix
    [V, D] = eig(covmat);
    D = diag(D);
    %the standard deviation of the posterior in model space is the square
    %root of each eigenvalue along its associated eigenvector.
    sig = D.^(1/2);
    medsig = median(sig);
    %initialize variables for modified eigenvalues
    Dmod = zeros(1,length(D));
    Dmodneg = zeros(1,length(D));
    %the modified eigenvalues will be those at which the likelihood is ~1/2
    %of the maximum likelihood, assuming a mulitvariate Gaussian
    %distribution. The logLkfn here returns the natural log of the
    %likelihood of the model * (-1).
    maxLk = logLkfn(opt);
    targetLk = maxLk + 1/2;
        
    for i = 1:length(sig)
        %step1--find two integer multiples of the eigenvalue which bracket
        %the target likelihood
        fa = 0;
        fb = 0;
        a = -1;
        b = 0;
        
        %set up system to check between 0 and 1 sigma, then after that move
        %to an exponential iteration system
        
        %boolean indicating whether to perform linear check between 0 and 1
        uselin = true;
        
        %check between 0 and 1
        perta = V(:,i) * medsig * 0;
        trya = opt + perta';
        pertb = V(:,i) * medsig * 1;
        tryb = opt + pertb';
            
        fa = logLkfn(trya)-targetLk;
        fb = logLkfn(tryb)-targetLk;
        
        %iterate exponentially if desired value not between 0 and 1
        while(fa*fb>=0)
            uselin = false;
            a = a + 1;
            b = b + 1;            
            
            perta = V(:,i) * medsig * exp(a);
            trya = opt + perta';
            pertb = V(:,i) * medsig * exp(b);
            tryb = opt + pertb';
            
            fa = logLkfn(trya)-targetLk;
            fb = logLkfn(tryb)-targetLk;
            
        end
        
        %step2--bisection algorithm
        
        %case of desired value lying between 0 and 1.
        if uselin
            a = 0;
            b = 1;
            err = 1;
            p = (a+b)/2;
            while(err>1e-7)
                p = (a+b)/2;
                pertp = V(:,i) * medsig * p;
                tryp = opt + pertp';
                fp = logLkfn(tryp)-targetLk;
                err = abs(fp);
                if(fa*fp<0)
                    b=p;
                else
                    a=p;
                end
            end
            Dmod(i) = (p*medsig)^2;
        %cases where desired value is >1 lie in exponential realm.
        else
            err = 1;
            p = (a+b)/2;
            while(err>1e-7)
                p = (a+b)/2;
                pertp = V(:,i) * medsig * exp(p);
                tryp = opt + pertp';
                fp = logLkfn(tryp)-targetLk;
                err = abs(fp);
                if(fa*fp<0)
                    b=p;
                else
                    a=p;
                end
            end
            Dmod(i) = (exp(p)*medsig)^2;
        end
    end
    
    %repeat above, but perturb in negative direction
    
    for i = 1:length(sig)
        %step1--find two integer multiples of the eigenvalue which bracket
        %the target likelihood
        fa = 0;
        fb = 0;
        a = -1;
        b = 0;
        
        %set up system to check between 0 and 1 sigma, then after that move
        %to an exponential iteration system
        
        %boolean indicating whether to perform linear check between 0 and 1
        uselin = true;
        
        %check between 0 and 1
        perta = V(:,i) * medsig * 0;
        trya = opt + perta';
        pertb = V(:,i) * medsig * -1;
        tryb = opt + pertb';
            
        fa = logLkfn(trya)-targetLk;
        fb = logLkfn(tryb)-targetLk;
        
%        keyboard;
        
        %iterate exponentially if desired value not between 0 and 1
        while(fa*fb>=0)
            uselin = false;
            a = a + 1;
            b = b + 1;            
            
            perta = V(:,i) * medsig * -exp(a);
            trya = opt + perta';
            pertb = V(:,i) * medsig * -exp(b);
            tryb = opt + pertb';
            
            fa = logLkfn(trya)-targetLk;
            fb = logLkfn(tryb)-targetLk;
            
            %-b
        end
        
        %step2--bisection algorithm
        %case of desired value lying between 0 and 1.
        if uselin
            a = 0;
            b = -1;
            err = 1;
            p = (a+b)/2;
            while(err>1e-7)
                p = (a+b)/2;
                pertp = V(:,i) * medsig * p;
                tryp = opt + pertp';
                fp = logLkfn(tryp)-targetLk;
                err = abs(fp);
                if(fa*fp<0)
                    b=p;
                else
                    a=p;
                end
            end
            Dmod(i) = (p*medsig)^2;
        %cases where desired value is >1 lie in exponential realm.
        else
            err = 1;
            p = (a+b)/2;
            while(err>1e-7)
                p = (a+b)/2;
                pertp = V(:,i) * medsig * -exp(p);
                tryp = opt + pertp';
                fp = logLkfn(tryp)-targetLk;
                err = abs(fp);
                if(fa*fp<0)
                    b=p;
                else
                    a=p;
                end
            end
            Dmod(i) = (exp(p)*medsig)^2;
        end
    end
    
    %set Dcheck to the mean distance from the maximum likelihood model
    %along each eigenvector at which the likelihood decreases to half the
    %maximum likelihood (there are two directions, positive and negative,
    %along each eigenvector).
    Dcheck = (Dmod+Dmodneg)/2;
    %solve for the matrix that has the same basis as the original
    %covariance matrix, but with the modified eigenvalues.
    B = V*diag((Dmod+Dmodneg)/2);
    y = (V'\B')';
catch
    keyboard;
end
end