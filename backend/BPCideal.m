function y = BPCideal(x1, x2, f1, f2)
%function calculates expected BPC value according to Eqn S7 in the text--x1 = N1/(N1+N2) where N1 and N2 are the sample sizes of samples 1 and 2.  f1 and f2 are the shared proportions of the compared populations.
    if(f1==0||f2==0||x1==0||x2==0)
        y = 0;
    else
        lLQ = x1*(1-f1)*log(x1) + (x1*f1+x2*f2)*log(x1*f1+x2*f2) + x2*(1-f2)*log(x2) - x1*f1*log(f1) - x2*f2*log(f2);
        lLQdiff = x1*log(x1) + x2*log(x2);
        y = 1-(lLQ/lLQdiff);
    end
end
