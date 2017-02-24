function y = evalBPC(folderpath, varargin)

%function evaluates the Bayesian Population Correlation (BPC) between all
%possible sample pairs of a set of zircon ages in folderpath/ whose
%posterior samples (Markov chains) have already been generated. The BPC
%values are displayed graphically.

%PARAMETERS
%folderpath - path of the directory where the zircon sample datafiles (age
%& uncertainty) are stored.  Markov chains should have already been
%generated using makechains() targeting this same path.

%varargin{1} - array of indices for desired sample order for BPC display
%table. Entering [1 2 3 4 5 ...] results in display in default order
%(alphabetical as retrieved by the dir() MATLAB function). Re-ordering
%indices changes the display order of the samples.

%OUTPUT
%y is an NxNx2 matrix, where N is the number of zircon samples or indices
%input into varargin{1}. y(:,:,1) contains the BPC values between each
%sample pair, as displayed in the figure created by this function, and
%y(:,:,2) contains the 1 sigma BPC uncertainty (standard deviation of
%log likelihood values of the Markov chain accepted models propagated
%through formula for BPC).

    %set flat variable for user-specified sample order.
    close all;
    if size(varargin,1)>0
        ORDERSPEC=1;
    else
        ORDERSPEC=0;
    end

    %retrieve sample names and initialize variables.
    fnames = dir(strcat(folderpath,'*.csv'));
    numfids = length(fnames);
    
    sLk = zeros(1,numfids);
    sstd = zeros(1,numfids);
    jLk = zeros(numfids);
    jstd = zeros(numfids);
    dLk = zeros(numfids);
    dstd = zeros(numfids);
    ssize = zeros(1,numfids);
    BPCm = ones(numfids);
    BPCs = zeros(numfids);
    
    %retrieve sample sizes from data files w/ age & uncertainties.
    for i = 1:numfids
        file = importdata(strcat(folderpath,fnames(i).name));
        if isstruct(file)
            ssize(i) = size(file.data,1);
        else
            ssize(i) = size(file,1);
        end
    end
    
    %read the mean and standard deviation of log likelihood of Markov chain
    %models of each individual sample 
    for i = 1:numfids
        logLk = csvread(strcat(folderpath,'chains/',fnames(i).name(1:end-4),'logLk.csv'),0);
        sLk(i) = mean(logLk);
        sstd(i) = std(logLk);
    end

    %read the mean and standard deviation of log likelihood of Markov chain
    %models of the 'joint' samples constructed from each possible
    %inter-sample combination and calculate the differences in likelihood
    %between the joint samples and the two combined samples on their own.
    for i = 1:numfids
        for j = i+1:numfids
            logLk = csvread(strcat(folderpath,'chains/',fnames(i).name(1:end-4),'_',fnames(j).name(1:end-4),'logLk.csv'),0);
            jLk(i,j) = mean(logLk);
            jstd(i,j) = std(logLk);
            dLk(i,j) = jLk(i,j) - sLk(i) - sLk(j);
            dstd(i,j) = sqrt(jstd(i,j)^2 + sstd(i)^2 + sstd(j)^2);
        end
    end
    
    %calculate BPC values
    for i = 1:numfids
        for j = i+1:numfids
            ssize1 = ssize(i);
            ssize2 = ssize(j);
            ssizetot = ssize1 + ssize2;
            BPCm(i,j) = 1 - dLk(i,j)/(ssize1*log(ssize1/ssizetot) + ssize2*log(ssize2/ssizetot));
            BPCm(j,i) = BPCm(i,j);
            BPCs(i,j) = dstd(i,j)/(-1*(ssize1*log(ssize1/ssizetot) + ssize2*log(ssize2/ssizetot)));
            BPCs(j,i) = BPCs(i,j);
        end
    end
    
    %reorder samples based on user input or based on which samples show the
    %greatest to least total BPC values combined with all other samples.
    sampnames = cell(size(BPCm,1),1);

    BPCsumsq = sum(BPCm.^2,2);
    if ~ORDERSPEC
        [s, idx] = sort(BPCsumsq);
    else
        idx = varargin{1,1};
    end

    for i = 1:length(idx)
        sampnames{i} = fnames(idx(i)).name;
    end

    BPCm = BPCm(idx,:);
    BPCm = BPCm(:,idx);

    BPCs = BPCs(idx,:);
    BPCs = BPCs(:,idx);

    %display BPC data in figure
    figure(1);
    image(BPCm.*64);

    set(gca, 'YTick', 1:size(BPCm,1), 'YTickLabel', sampnames);
    set(gca, 'XTick', 1:size(BPCm,1), 'XTickLabel', sampnames);

    for i = 1:length(idx)
        for j = 1:length(idx)
            if i==j
                text(j,i,'-','FontSize',12,'FontWeight','bold');
            else
                l_str = num2str(round(BPCm(i,j)*100)/100);
                l_err = num2str(round(BPCs(i,j)*100)/100);
                cut = length(l_str);
                cuterr = length(l_err);
                if cut>4
                    cut = 4;
                end
                if cuterr>4
                    cuterr = 4;
                end
                l_str = l_str(1:cut);
                l_err = l_err(2:cuterr);
                text(j-0.26*(cut/4),i-0.15,l_str,'FontSize',12,'FontWeight','bold');
                text(j-0.26*(cuterr/4),i+0.2,strcat('±',l_err),'FontSize',10,'FontWeight','normal');
            end
        end
    end
    
    %set y to matrix of BPC values and uncertainties.
    y = zeros(size(BPCm,1), size(BPCm,1), 2);
    y(:,:,1) = BPCm;
    y(:,:,2) = BPCs;
end