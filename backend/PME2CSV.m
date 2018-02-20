function z = PME2CSV(folderpath,name,x_min,x_max,xvals,logPDFs,savepath,savename,varargin)
    
%OUTPUT
%z is a dummy variable

%PARAMETERS
%folderpath--the path where the .csv file of the zircon ages is located.
%Same as the folderpath input for makePME() and evalBPC().

%name--the name of the sample for which display is desired.  This is the
%same as the name of the sample's .csv file, except with no extension.

%x_min, x_max--these are the minimum and maximum x values used to generate
%the chains for the sample to be displayed (same as input into
%makePME()).

%xvals--the age values for which probability values will be written to CSV.

%logPDFs--boolean that determines whether PDFs will be shown on log p scale
%or with linear p scale.

%logoutput--boolean that determines whether the output CSV will be gridded
%in terms of linear or log age.

%varargin(1)--number of cores

%varargin(2)--max. number of PDFs to export

    if size(varargin,2)>0
        corespec = varargin{1};
    else
        corespec = 0;
    end
    
    if size(varargin,2)>1
        thin = varargin{2};
    else
        thin = 0;
    end
    
    %A corespect value of -1 indicates not to close and reopen the parallel
    %pool
    if corespec~=-1
        delete(gcp('nocreate'));
    end
    
    %start parpool
    if corespec == 0
        defaultProfile = parallel.defaultClusterProfile;
        myCluster = parcluster(defaultProfile);
        parpool(myCluster);
    elseif corespec~=-1
        defaultProfile = parallel.defaultClusterProfile;
        myCluster = parcluster(defaultProfile);
        parpool(myCluster,corespec);
    end    

    %log transform x_min and x_max and xvals
    x_min = log(x_min);
    x_max = log(x_max);
    xvalslin = xvals;
    xvals = log(xvals);
    
    %set knot locations for spline
    N_coefs=50; %by default
    knots = [x_min x_min x_min:(x_max-x_min)/(N_coefs-2):x_max x_max x_max];
    %read pchain
    pchain = csvread(strcat(folderpath,'chains/',name,'chain.csv'),0);
    
    if thin>0
        idx = randperm(size(pchain,1),thin);
        pchain = pchain(idx,:);
    end
    
    %set x values at which to evaluate each PDF, dictated by x_min,
    %x_max, and x_res

    PDFs = zeros(size(pchain,1),length(xvals));
    
    %generate logPDF or linear PDF function values at indicated x-values 
    %for all probability chains

    h = parfor_progressbar(size(PDFs,1),'Evaluating PDFs...');
    parfor i = 1:size(PDFs,1)
        if ~logPDFs
            PDFs(i,:) = exp(fnval(spmak(knots,pchain(i,:)),xvals));
        else
            PDFs(i,:) = fnval(spmak(knots,pchain(i,:)),xvals);
        end
        h.iterate();
    end
    close(h);
    
    %renormalize PDF area to match the x scale
    if ~logPDFs
        PDFs = PDFs * (x_max-x_min)/(exp(x_max)-exp(x_min));
    else
        PDFs = PDFs + log((x_max-x_min)/(exp(x_max)-exp(x_min)));
    end
    
    csvwrite(strcat(savepath,savename),[xvalslin';PDFs]');
    
    z=0;
end