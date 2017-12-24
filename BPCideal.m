function y = BPCideal(x1, x2, f1, f2)
    if(f1==0||f2==0||x1==0||x2==0)
        y = 0;
    else
        lLQ = x1*(1-f1)*log(x1) + (x1*f1+x2*f2)*log(x1*f1+x2*f2) + x2*(1-f2)*log(x2) - x1*f1*log(f1) - x2*f2*log(f2);
        lLQdiff = x1*log(x1) + x2*log(x2);
        y = 1-(lLQ/lLQdiff);
    end
end