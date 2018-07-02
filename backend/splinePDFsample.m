function y = splinePDFsample(knots,coefs,n,varargin)
    %PARAMETERS
    %knots are the knots of a spline function that is to be interpreted as
    %a log-scale PDF
    %coefs are the spline coefficients for a log-scale PDF
    
    %n is the number of desired samples from the PDF.
    
    %varargin{1} is optionally the number of sets of samples to obtain (n
    %is the number of ages, N is the number of sets of ages or age samples)
    
    %OUTPUT
    %function draws samples randomly from the range of the spline function,
    %and for each one generates a random number from [0,max(PDFvalue)]. If
    %this second random number falls beneath the PDF curve at the x value
    %of the first number, then the first is accepted into the sample.  In
    %this way, the PDF is sampled accurately.
    
    if size(varargin,2)>0
        N = varargin{1};
    else
        N = 1;
    end

    try
    sample = zeros(n,N);
    PDF = spmak(knots,coefs);
    PDFmax = max(exp(fnval(PDF,[min(knots):range(knots)/(10*length(knots)):max(knots)])));
    
    for k = 1:N
        curidx = 1;
        while curidx<n+1
            x = unifrnd(min(knots),max(knots),n,1);
            y = unifrnd(0,2*PDFmax,n,1);
            %2 * PDFmax for wiggle room (in case the true maximum isn't caught
            %by the discretization used to calculate PDFmax above

            toadd = find(log(y)<fnval(PDF,x));
            if(curidx+length(toadd)>n+1)
                toadd = toadd(1:n+1-curidx);
            end

            if(length(toadd)>0)
                sample(curidx:curidx+length(toadd)-1,k) = x(toadd);
                curidx = curidx + length(toadd);
            end
        end
        
    end

    y = sample;
    catch e
        fprintf(1,'The identifier was:\n%s',e.identifier);
        fprintf(1,'There was an error! The message was:\n%s',e.message);
        keyboard
    end
end