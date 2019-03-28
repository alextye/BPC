function canceled = makePME(folderpath, x_min, x_max, varargin)

%function generates representative samples of the posteriors of the set of
%zircon samples contained within folderpath, as well as the 'joint' samples
%that represent all possible inter-sample combinations of the samples in
%folderpath

%OUTPUT
%y is a dummy variable. The sample of the posterior, found via Markov 
%chain, are all saved in the directory 'folderpath/chains/', with the 
%naming convention 'Sample_namechain.csv' for the spline coefficients of 
%each model and 'Sample_namelogLk.csv' for the natural log likelihood
%values of each model. The Markov chain spline coefficients and natural log
%likelihood values of the inter-sample combinations are saved in the same
%folder with naming conventions 'Sample_name1Sample_name2chain.csv' and
%'Sample_name1Sample_name2logLk.csv'. 'Sample_name', 'Sample_name1', and
%'Sample_name2' are taken from the names of the files in folderpath.

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

%varargin{1} is optionally the number of cores to use for processing
%varargin{2} is optionally the flag that tells the routine whether to
%account for analytical uncertainties on grain ages.
%varargin{3} optionally flags the non-Bayesian 'quick mode'
%varargin{4} flags whether to overwrite existing PME files (default is not)
%varargin{5} is optionally the number of models to retain from the inferred
    %PMEs

    delete(gcp('nocreate'));
    
    %ensure the parallel pool is closed upon completion
    cancelFutures = onCleanup(@() delete(gcp('nocreate')));
    canceled = 0;

    if size(varargin,2)>4
        maxmodelspec = varargin{5};
    else
        maxmodelspec = 0;
    end
    
    if size(varargin,2)>3
        OWFLAG = varargin{4};
    else
        OWFLAG = 0;
    end
    
    if size(varargin,2)>2
        QUICKMFLAG = varargin{3};
    else
        QUICKMFLAG = 0;
    end
    
    if size(varargin,2)>1
        AUFLAG = varargin{2};
    else
        AUFLAG = 1;
    end
    
    if size(varargin,2)>0
        corespec = varargin{1};
    else
        corespec = 0;
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
        if ~AUFLAG
            lsampages{i}(:,2) = 0;
        end
    end
    
    
    N_coefs = 100;
    %N_coefs = 250;
    p_sig = 13.5;
    lxmin = log(x_min);
    lxmax = log(x_max);
    N_pts = 1000;
    delta = N_coefs/10;
    knotbreak = floor(5/7*(N_coefs-2));
    %knots = [lxmin lxmin lxmin:(lxmax-lxmin)/(N_coefs-2):lxmax lxmax lxmax];
    
    %construct knots such that they are log distributed <1000 Ma and linear
    %distributed >1000 Ma, mirroring the distribution of zircon U-Pb age 
    %errors.
    knotslog = [lxmin lxmin lxmin:(log(1000)-lxmin)/knotbreak:log(1000)];
    knotslin = [1000:(x_max-1000)/(N_coefs-2-knotbreak):x_max x_max x_max];
    knots = [knotslog log(knotslin(2:end))];
    
    %create directory for the Markov chain results
    out = evalc('mkdir(strcat(folderpath,''chains/''))');
    %create directory for log files
    out = evalc('mkdir(strcat(folderpath,''log/''))');

    %set up index "guide" for better parallelization below
    idxguide = zeros(numfids);
    k = numfids;
    for i = 1:numfids
        for j = i+1:numfids
            k = k+1;
            idxguide(i,j) = k;
        end
    end

    h = waitbar(0,'Inferring PMEs... (can take minutes to hours)','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

    %iterate through the zircon sample datafiles in folderpath, generating
    %a Markov chain (posterior sample of possible probability models) for
    %each zircon sample.
    
    numActive = 0;
    isactive = zeros(1,k);

    for l = 1:k
                
        if l<=numfids
            i = l;
            %test for whether the Markov chain already exists; this allows the
            %script to be terminated and restarted with minimal consequence.
            %if the user has chosen to overwrite existing files, then
            %perform the analysis no matter whether the file already exists
            %or not.
            nametest = dir(strcat(folderpath,'chains/',fnames(i).name(1:end-4),'logLk.csv'));
            if(size(nametest,1)==0|OWFLAG)
                
                %start parpool if not already open
                if isempty(gcp('nocreate'))
                    if corespec == 0
                        defaultProfile = parallel.defaultClusterProfile;
                        myCluster = parcluster(defaultProfile);
                        parpool(myCluster);
                    else
                        defaultProfile = parallel.defaultClusterProfile;
                        myCluster = parcluster(defaultProfile);
                        parpool(myCluster,corespec);
                    end
                end
                
                numActive = numActive + 1;
                isactive(l) = 1;
                filename = strcat(folderpath,'log/',fnames(i).name(1:end-4),'log.txt');
                
                %add the operations to the parallelized function queue.
                F(numActive) = parfeval(@chain_const,3,N_coefs,knots,lsampages{i},lxmin,lxmax,p_sig,delta,N_pts,filename,QUICKMFLAG);
            end
        else

    %iterate through the zircon sample datafiles in folderpath, generating
    %a Markov chain (posterior sample of possible probability models) for
    %the 'joint' sample of each inter-sample combination.

            [i, j] = find(idxguide==l);
            nametest = dir(strcat(folderpath,'chains/',fnames(i).name(1:end-4),'_',fnames(j).name(1:end-4),'logPost.csv'));
            if(size(nametest,1)==0|OWFLAG)
                numActive = numActive + 1;
                isactive(l) = 1;
                filename = strcat(folderpath,'log/',fnames(i).name(1:end-4),'_',fnames(j).name(1:end-4),'log.txt');
                
                %add the operations to the parallelized function queue.
                F(numActive) = parfeval(@chain_const_joint,3,N_coefs,knots,lsampages{i},lsampages{j},lxmin,lxmax,p_sig,delta,N_pts,filename,QUICKMFLAG);
            end
        end
        
    end
    
    waitbar(1-numActive/k);
    
    activeIdx = find(isactive);
    l = 0;

    %retrieve the results from the parallelized function queue while also
    %continuing to check for whether the user has canceled the operation.
    while l<numActive
        
        if getappdata(h,'canceling')
            canceled = 1;
            break
        end
        
        [Fidx, pchain, logLk, logPost] = fetchNext(F,2);
        
        %thin chain for writing to disk if specified
        if length(logLk) > maxmodelspec & maxmodelspec > 0
            %preserve max. likelihood model at head of each chain
            pchain1 = pchain(1,:);
            logLk1 = logLk(1);
            logPost1 = logPost(1);
            
            writeidx = floor([1:maxmodelspec] * length(logLk)/maxmodelspec);
            pchain = pchain(writeidx,:);
            logLk = logLk(writeidx);
            logPost = logPost(writeidx);
            
            pchain(1,:) = pchain1;
            logLk(1) = logLk1;
            logPost(1) = logPost1;
        end
        
        if ~isempty(Fidx)
            completedIdx = activeIdx(Fidx);
        
            if completedIdx <= numfids
                i = completedIdx;
                dlmwrite(strcat(folderpath,'chains/',fnames(i).name(1:end-4),'chain.csv'),knots);
                dlmwrite(strcat(folderpath,'chains/',fnames(i).name(1:end-4),'chain.csv'),pchain,'-append');
                csvwrite(strcat(folderpath,'chains/',fnames(i).name(1:end-4),'logLk.csv'),logLk);
                csvwrite(strcat(folderpath,'chains/',fnames(i).name(1:end-4),'logPost.csv'),logPost);
            else
                [i, j] = find(idxguide==completedIdx);
                dlmwrite(strcat(folderpath,'chains/',fnames(i).name(1:end-4),'_',fnames(j).name(1:end-4),'chain.csv'),knots);
                dlmwrite(strcat(folderpath,'chains/',fnames(i).name(1:end-4),'_',fnames(j).name(1:end-4),'chain.csv'),pchain,'-append');
                csvwrite(strcat(folderpath,'chains/',fnames(i).name(1:end-4),'_',fnames(j).name(1:end-4),'logLk.csv'),logLk);
                csvwrite(strcat(folderpath,'chains/',fnames(i).name(1:end-4),'_',fnames(j).name(1:end-4),'logPost.csv'),logPost);
            end
            l = l + 1;
            waitbar(1-(numActive-l)/k);
        end
        
    end
    
    delete(h);
    delete(gcp('nocreate'));
end