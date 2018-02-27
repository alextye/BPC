function y = BPCunc(folderpath, x_min, x_max, varargin)

%in order to estimate BPC uncertainty, function generates random subsamples
%from PMEs of samples contained in folderpath. These random subsamples are
%generated by selecting PDFs at random from the PMEs, and then randomly
%sampling each of those PDFs. The result is a set of synthetic zircon age
%samples each with a number of ages equal to the size of the true detrital
%zircon age sample for which the PME was inferred.  The synthetic age
%samples are then compared to one another using Maximum Likelihood
%Correlation, as discussed in the text, and spread in the results is used
%as an estimate of BPC uncertainty.

%OUTPUT
%y is a dummy variable. The script creates datafiles saved in the directory
%folderpath/unc/. These files contain the synthetic age samples discussed
%above as well as the likelihood of the maximum likelihood probability
%model for each synthetic sample.

%PARAMETERS
%folderpath is the path of the directory where the zircon age data of the
%zircon samples is stored. Each zircon sample should have one .csv file,
%named after the sample name, and with the .csv extension:
%'samplename.csv'. Each .csv file should follow the format below:
%
%age1,unc1,
%age2,unc2,
%age3,unc3,
%...
%where age1, age2, and age3 correspond to measured ages of grains 1, 2, and
%3, and unc1, unc2, and unc3 correspond to analytical uncertainties of
%grains 1, 2, and 3.

%x_min and x_max are the minimum and maximum desired x values for the
%modeling. Decreasing the range of x_min and x_max results in greater
%splining resolution over the specified range, though any grain ages lower
%than x_min or greater than x_max are not considered in the analysis. NOTE: 
%We suggest that BPC not be evaluated between posterior samples generated 
%with differing x_min and x_max values, so x_min and x_max should be as
%inclusive as necessary for all samples to be analyzed using BPC. NOTE 2:
%Our splining procedure uses a natural logarithmic age scale, so x_min
%values of 0 are prohibited and values of <1 are not suggested.

%varargin is optionally the number of cores to use for processing

    delete(gcp('nocreate'));

    if size(varargin,2)>0
        corespec = varargin{1};
    else
        corespec = 0;
    end

    %start parpool
    if corespec == 0
        defaultProfile = parallel.defaultClusterProfile;
        myCluster = parcluster(defaultProfile);
        parpool(myCluster);
    else
        defaultProfile = parallel.defaultClusterProfile;
        myCluster = parcluster(defaultProfile);
        parpool(myCluster,corespec);
    end

    %retrieve the names of zircon sample data files.
    fnames = dir(strcat(folderpath,'*.csv'));
    
    numfids = length(fnames);
    sampages = cell(1,numfids);
    lsampages = cell(1,numfids);
    
    %retrieve zircon sample data
    for i = 1:numfids
        file = importdata(strcat(folderpath,fnames(i).name));
        
        if isstruct(file)
            sampages{i} = file.data;
        else
            sampages{i} = file;
        end
        
        
        col = size(sampages{i},2);
        if(col>2)
            sampages{i} = [sampages{i}(:,col-1) sampages{i}(:,col)];
        end
        lsampages{i} = sampages{i};
        lsampages{i}(:,1) = log(sampages{i}(:,1));
        lsampages{i}(:,2) = sampages{i}(:,2)./sampages{i}(:,1);
    end
    
    %set parameters of resampling for uncertainty calculation--M is the
    %number of models selected for each PME, N is the number of samples
    %selected from each model
    M = 5;
    N = 20;
    
    %set parameters for splining and MCMC. See comments of chain_const.m
    %for discussion.
    N_coefs = 50;
    p_sig = 13.5;
    lxmin = log(x_min);
    lxmax = log(x_max);
    N_pts = 1000;
    delta = N_coefs/10;
    knots = [lxmin lxmin lxmin:(lxmax-lxmin)/(N_coefs-2):lxmax lxmax lxmax];
    
    %create working directory for operations used to estimate uncertainty
    evalc('mkdir(strcat(folderpath,''unc/''))');
    evalc('mkdir(strcat(folderpath,''unc/log/''))');
    
    %set up index "guide" for better parallelization below
    idxguide = zeros(numfids);
    k = 0;
    for i = 1:numfids
        for j = i+1:numfids
            k = k+1;
            idxguide(i,j) = k;
        end
    end

    h = parfor_progressbar(k+numfids,'Estimating uncertainties... (can take minutes to hours)');
    %iterate through the zircon sample datafiles in folderpath, generating
    %synthetic samples for each real zircon sample by resampling PMEs.
    
    parfor i = 1:numfids
        %test for whether the resamples already exist; if so, bypass this
        %iteration. this allows the
        %script to be terminated and restarted with minimal consequence.
        nametest = dir(strcat(folderpath,'unc/',fnames(i).name(1:end-4),'resamples.csv'));
        
        if(size(nametest,1)==0)
        
            fileID = fopen(strcat(folderpath,'unc/log/',fnames(i).name(1:end-4),'log.txt'),'w');
            pchain = csvread(strcat(folderpath,'chains/',fnames(i).name(1:end-4),'chain.csv'));
            idx = ceil(size(pchain,1)*rand(M,1));
            n = size(lsampages{i},1);
            samples = zeros(n,N*M);

            %generate N random samples of n grains from each of the M
            %models drawn from each PME.
            for m = 1:M
                samples(:,[(m-1)*N+1:m*N]) = splinePDFsample(knots,pchain(idx(m),:),n,N);
            end
            %calculate the maximum likelihood model and logLk for each
            %sample
            logLk = zeros(1,N*M);
            for m = 1:M*N
                [pchain, logLk(m)] = chain_const_noMCMC(N_coefs, [samples(:,m) zeros(n,1)], lxmin, lxmax, p_sig, delta, N_pts, fileID);
            end
            
            %write the output in the following format: n+1 x (M*N) 
            %matrix, each column one random subsample from the PME,
            %each line a grain age, except first line logLk of the
            %maximum likelihood model for each random subsample
            csvwrite(strcat(folderpath,'unc/',fnames(i).name(1:end-4),'resamples.csv'),[logLk; samples]);
            
            fclose('all');
        end
        h.iterate();
    end

    %iterate through the possible pairs of zircon samples in folderpath,
    %comparing the synthetic samples for each pair.
    parfor l = 1:k
        [i, j] = find(idxguide==l);
        nametest = dir(strcat(folderpath,'unc/',fnames(i).name(1:end-4),'_',fnames(j).name(1:end-4),'jointML.csv'));
        if(size(nametest,1)==0)
            %read the files corresponding to the PME subsamples for
            %each detrital zircon sample corresponding to i, j
            resamples1 = csvread(strcat(folderpath,'unc/',fnames(i).name(1:end-4),'resamples.csv'),1);
            resamples2 = csvread(strcat(folderpath,'unc/',fnames(j).name(1:end-4),'resamples.csv'),1);
            fileID = fopen(strcat(folderpath,'unc/log/',fnames(i).name(1:end-4),'_',fnames(j).name(1:end-4),'log.txt'),'w');
            
            n1 = size(lsampages{i},1);
            n2 = size(lsampages{j},1);

            %initialize logLk array for the maximum likelihood models
            %of the joint samples
            jlogLk = zeros(N*M,1);

            %find the maximum likelihood model and maximum logLk for
            %joint samples composed of pairs of the PME random
            %subsamples generated above
            for m = 1:N*M
                
                [jointchain, jlogLk(m)] = chain_const_noMCMC(N_coefs,[resamples1(:,m) zeros(n1,1);resamples2(:,m) zeros(n2,1)],lxmin,lxmax,p_sig,delta,N_pts,fileID);
                
            end
            
            %write the maximum logLk values for the pair of
            %detrital zircon samples out to a file.
                
            csvwrite(strcat(folderpath,'unc/',fnames(i).name(1:end-4),'_',fnames(j).name(1:end-4),'jointML.csv'),jlogLk);
            
            fclose('all');
        end
        h.iterate();
    end
    close(h);
    delete(gcp('nocreate'));
%    h = msgbox('BPCunc complete.');
end