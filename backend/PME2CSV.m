function z = PME2CSV(folderpath,name,x_min,x_max,Npts,logage,logPDFs,savepath,savename,varargin)
    
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
    
    %A corespec value of -1 indicates not to close and reopen the parallel
    %pool, for instance if multiple samples are being queried.
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
    
    %set knot locations for spline
    %N_coefs=50; %by default
    %knots = [x_min x_min x_min:(x_max-x_min)/(N_coefs-2):x_max x_max x_max];
    %read pchain
    data = csvread(strcat(folderpath,'chains/',name,'chain.csv'),0);
    knots = data(1,:);
    %read pchain
    pchain = data([2:end],[1:end-3]);
    
    if thin>0 & thin<size(pchain,1)
        idx = randperm(size(pchain,1),thin);
        pchain = pchain(idx,:);
    end
        
    
    PDFs = zeros(size(pchain,1),Npts);
    integrals = zeros(size(pchain,1),1);
    
    
    %set x values at which to evaluate each PDF, dictated by x_min,
    %x_max, and x_res

    %set delx variable used to estimate PDF integrals
    xvals = zeros(1,Npts);
    if logage
        delx = (x_max-x_min)/Npts;
        xvals = linspace(x_min+delx/2,x_max-delx/2,Npts);
    else
        delx = (exp(x_max)-exp(x_min))/Npts;
        xvals = log(linspace(exp(x_min)+delx/2,exp(x_max)-delx/2,Npts));        
    end
    
    %keyboard
    
    %generate logPDF or linear PDF function values at indicated x-values 
    %for all probability chains
    
    h = parfor_progressbar(size(PDFs,1),'Evaluating PDFs...');
    parfor i = 1:size(PDFs,1)
        if ~logPDFs & ~logage
            %if linear age is requested, divide the PDF values (valid for
            %log age scale) by the exponentiated x values, which results in
            %properly scaled peak heights
            PDFs(i,:) = fnval(spmak(knots,pchain(i,:)),xvals);
            PDFs(i,:) = exp(PDFs(i,:)-xvals);
        elseif ~logPDFs & logage
            PDFs(i,:) = exp(fnval(spmak(knots,pchain(i,:)),xvals));         
            %if linear age is requested, divide the PDF values (valid for
            %log age scale) by the exponentiated x values, which results in
            %properly scaled peak heights
%             if ~logage
%                 PDFs(i,:) = PDFs(i,:)./exp(xvals);
%             end
            
        else
            PDFs(i,:) = fnval(spmak(knots,pchain(i,:)),xvals);
            if ~logage
                PDFs(i,:) = PDFs(i,:)-xvals;
            end
        end
        
        h.iterate();
    end
    close(h);
    
	%if PDFs in linear age space are requested, exponentiate the x values for output.
    if ~logage
    	xvals = exp(xvals);
    end
    
    %calculate Riemann summation approximate integrals
    if ~logPDFs
        %integrals = delx * (sum(PDFs(:,[2:Npts-1]),2) + .5 * (PDFs(:,1)+PDFs(:,Npts)));
        integrals = delx * sum(PDFs,2);
    else
        %integrals = delx * (sum(exp(PDFs(:,[2:Npts-1])),2) + .5 * (exp(PDFs(:,1))+exp(PDFs(:,Npts))));
        integrals = delx * sum(exp(PDFs),2);
    end
    
%    keyboard
    
    %write header to file
    fid = fopen(strcat(savepath,savename),'w'); 
    fprintf(fid,'%s\n','Upper content line shows integrals (Riemann summation approximation method) of each PDF. Below this line the first column is x values and the successive columns are PDF values at the given x values. Integral values in first line correspond to PDF values directly below.  If integrals deviate from 1 too much increase the ''number of points'' field in the GUI.');
    fprintf(fid,'%s\n','Integrals');
    fclose(fid);
    %write data to end of file
    dlmwrite(strcat(savepath,savename),[0 integrals'],'-append');
    fid = fopen(strcat(savepath,savename),'a'); 
    fprintf(fid,'%s\n','x values, PDF values');
    fclose(fid);
    dlmwrite(strcat(savepath,savename),[xvals;PDFs]','-append');

%    csvwrite(strcat(savepath,savename),[xvals;PDFs]');
    
    z=0;
end