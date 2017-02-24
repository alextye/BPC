%ALEX TYE
%22 NOV 2015
function y = BPC2frac(BPCm, BPCs, N1, N2, varargin)
%function BPC2frac takes in BPC values and related arguments and
%determines possible pairs of (f1, f2; the overlapping proportions of
%two detrital zircon populations) that could have produced the observed
%BPC value.
    
%PARAMETERS
%BPCm is the BPC value
%BPCs is the BPC uncertainty (standard deviation of
%log likelihood values of the Markov chain accepted models propagated
%through formula for BPC).
%N1, N2 are sample sizes of samples 1, 2

%varargin takes up to two additional options:
%varargin{1} argument is called SUPPLOT, which can be 0 or
%1. If 1, then plots associated with this function are suppressed.
        
%varargin{2} is f1. If f1 is specified, then the
%function finds values of f2 that could have produced the observed
%BPC values, and plots a PDF of potential values of f2. If f1 is
%specified, the function returns the [mean, standard_deviation] of
%the calculated f2 density function.

%OUTPUT
%function returns the mean and standard deviation of possible values of
%f2 if f1 is specified (see below), or returns [0 0] in the case of no
%f1 specification (the solution for (f1, f2) is complicated).
    
        
    %initiatlize variables and read in optional arguments
    x1 = N1/(N1+N2);
    x2 = N2/(N1+N2);
    f1SPEC = 0;
    SUPPLOT = 0;
    if size(varargin,2)>1
        f1 = varargin{1,2};
        SUPPLOT = varargin{1,1};
        f1SPEC = 1;
    elseif size(varargin,2)>0
        SUPPLOT = varargin{1,1};
    end
    
    if f1SPEC
        x = [0.01:0.001:0.99];
        z = zeros(size(x));
        inputdens = zeros(size(x));
        %calculate idealized BPC values for a range of possible f2 values
        %(stored in array x) and evaluate the probability density of
        %observing such BPC values given the BPCm and BPCs input as
        %function arguments.
        for i = 1:size(x,2)
            z(i) = BPCideal(x1,x2,f1,x(i));
            inputdens(i) = normpdf(z(i),BPCm,BPCs);
        end
        %plot the probability density of tested f2 values
        if(~SUPPLOT) plot(x,inputdens); end
        %find the center and standard deviation of the (assumed) normal
        %distribution of f2 by search
        mi = find(inputdens==max(max(inputdens)));
        si = find(abs(inputdens-normpdf(BPCm-BPCs,BPCm,BPCs))==min(abs(inputdens-normpdf(BPCm-BPCs,BPCm,BPCs))));
        y = [x(mi) abs(x(si)-x(mi))];
    else
        f1SPEC=0;
        %generate reference tables of f1,f2 values (called x,y) that will
        %have their BPC values calculated.
        x = [0.01:0.001:0.99];
        y = [-0.99:0.001:-0.01];
        ymap = [0.01:0.001:0.99];
        y = y';
        ymap = ymap';
        x = repmat(x,size(y,1),1);
        y = repmat(y,1,size(x,1));
        ymap = repmat(ymap,1,size(x,1));
        z = zeros(size(x));
        inputdens = zeros(size(x));
        %loop through reference tables calculating idealized BPC for given 
        %(f1,f2), and evaluating the probability density of observing such 
        %a value, given the BPCm and BPCs input as function arguments.
        for i = 1:size(x,1)
            for j = 1:size(y,2)
                z(i,j) = BPCideal(x1,x2,x(1,i),abs(ymap(j,1)));
                inputdens(i,j) = normpdf(z(i,j),BPCm,BPCs);
            end
        end
        if(~SUPPLOT)
            %plot the probability density of observing each idealized BPC
            %value given the BPCm and BPCs values input as function
            %arguments
            image(inputdens/max(max(inputdens)).*64);
            set(gca,'XTickLabel',{'0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9'});
            set(gca,'YTickLabel',{'0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9'});
        end
        y = [0 0];
    end

end