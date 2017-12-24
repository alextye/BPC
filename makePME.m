function y = makePME(folderpath, x_min, x_max, varargin)

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

%varargin is optionally the number of cores to use for processing

    delete(gcp('nocreate'));

    if size(varargin,2)>0
        corespec = varargin{1};
    else
        corespec = 0;
    end

    %start parpool
    if corespec == 0
        parpool;
    else
        parpool(corespec);
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
    
    
    %set parameters for splining and MCMC. See comments of chain_const.m
    %for discussion.
    N_coefs = 50;
    p_sig = 13.5;
    lxmin = log(x_min);
    lxmax = log(x_max);
    N_pts = 1000;
    delta = N_coefs/10;
    
    %create directory for the Markov chain results
    mkdir(strcat(folderpath,'chains/'));
    %create directory for log files
    mkdir(strcat(folderpath,'log/'));

    %set up index "guide" for better parallelization below
    idxguide = zeros(numfids);
    k = numfids;
    for i = 1:numfids
        for j = i+1:numfids
            k = k+1;
            idxguide(i,j) = k;
        end
    end

    h = parfor_progressbar(k,'Inferring PMEs... (can take minutes to hours)');
    %iterate through the zircon sample datafiles in folderpath, generating
    %a Markov chain (posterior sample of possible probability models) for
    %each zircon sample.
    parfor l = 1:k
        if l<=numfids
            i = l;
            %test for whether the Markov chain already exists; this allows the
            %script to be terminated and restarted with minimal consequence.
            nametest = dir(strcat(folderpath,'chains/',fnames(i).name(1:end-4),'logLk.csv'));
            if(size(nametest,1)==0)
                fileID = fopen(strcat(folderpath,'log/',fnames(i).name(1:end-4),'log.txt'),'w');
                [pchain, logLk] = chain_const(N_coefs,lsampages{i},lxmin,lxmax,p_sig,delta,N_pts,fileID);
                csvwrite(strcat(folderpath,'chains/',fnames(i).name(1:end-4),'chain.csv'),pchain);
                csvwrite(strcat(folderpath,'chains/',fnames(i).name(1:end-4),'logLk.csv'),logLk);
                fclose('all');
            end
        else

    %iterate through the zircon sample datafiles in folderpath, generating
    %a Markov chain (posterior sample of possible probability models) for
    %the 'joint' sample of each inter-sample combination.

            [i, j] = find(idxguide==l);
            nametest = dir(strcat(folderpath,'chains/',fnames(i).name(1:end-4),'_',fnames(j).name(1:end-4),'logLk.csv'));
            if(size(nametest,1)==0)
                fileID = fopen(strcat(folderpath,'log/',fnames(i).name(1:end-4),'_',fnames(j).name(1:end-4),'log.txt'),'w');
                [jointchain, jointlogLk] = chain_const_joint(N_coefs,lsampages{i},lsampages{j},lxmin,lxmax,p_sig,delta,N_pts,fileID);
                csvwrite(strcat(folderpath,'chains/',fnames(i).name(1:end-4),'_',fnames(j).name(1:end-4),'chain.csv'),jointchain);
                csvwrite(strcat(folderpath,'chains/',fnames(i).name(1:end-4),'_',fnames(j).name(1:end-4),'logLk.csv'),jointlogLk);
                fclose('all');
            end
        end
        h.iterate();
    end
    close(h);
    delete(gcp('nocreate'));
    h = msgbox('makePME complete.');
end