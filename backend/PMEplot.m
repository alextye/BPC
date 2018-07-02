function z = PMEplotlin(folderpath,name,x_min,x_max,x_res,y_res,fignum,showMLM,showdot,logPDFs,thin,varargin)
    
%OUTPUT
%z is a dummy variable

%PARAMETERS
%folderpath--the path where the .csv file of the zircon ages is located.
%Same as the folderpath input for makechains() and evalBPC().

%name--the name of the sample for which display is desired.  This is the
%same as the name of the sample's .csv file, except with no extension.

%x_min, x_max--these are the minimum and maximum x values used to generate
%the chains for the sample to be displayed (same as input into
%makechains()).

%x_res, y_res--these are the desired resolution of the plot for the x- and
%y-coordinates (higher number = finer resolution, longer runtime).

%fignum--this is the number of the figure to which to output the display

%showMLM--this is a boolean that determines whether to show the maximum
%likelihood probability model as a colored line in the resulting figure or
%not.

%showdot--boolean that determines whether a dotplot of the measured ages
%will be shown beneath the PME.

%logPDFs--boolean that determines whether PDFs will be shown on log p scale
%or with linear p scale.

%varargin(1)--number of cores

%varargin(2)--optional plot title.

%varargin(3)--indicates whether dot plot and maximum likelihood model
%should be plotted in new figures

    if size(varargin,2)>2
        showtitle=1;
        plottitle = varargin{3};
    else
        showtitle=0;
    end

    if size(varargin,2)>1
        newfig=varargin{2};
    else
        newfig=0;
    end
    
    if size(varargin,2)>0
        corespec = varargin{1};
    else
        corespec = 0;
    end
    
    %CORESPEC of -1 tell the script to use whatever parallel pool is
    %already open and not to close the pool at end of run.
    if(corespec~=-1)
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

    %log transform x_min and x_max
    x_min = log(x_min);
    x_max = log(x_max);

    %set knot locations for spline
    N_coefs=50; %by default
    knots = [x_min x_min x_min:(x_max-x_min)/(N_coefs-2):x_max x_max x_max];
    %read pchain
    pchain = csvread(strcat(folderpath,'chains/',name,'chain.csv'),0);
    
    %thin (select a subset of models) for plotting if user has specified to
    %do so.
    if thin > 0 & thin < size(pchain,1)
        %preserve MLM at head of pchain
        pchain1 = pchain(1,:);
        
        idx = randperm(size(pchain,1)-1,thin);
        pchain = pchain(idx,:);
        pchain(1,:) = pchain1;
    end
    
    %set x values at which to evaluate each PDF, dictated by x_min,
    %x_max, and x_res
    x = [x_min+(x_max-x_min)/(2*x_res):(x_max-x_min)/x_res:x_max-(x_max-x_min)/(2*x_res)];
    PDFs = zeros(size(pchain,1),x_res);
    
    %generate logPDF or linear PDF function values at indicated x-values 
    %for all probability chains

    h = parfor_progressbar(size(PDFs,1),'Evaluating PDFs...');
    parfor i = 1:size(PDFs,1)
        if ~logPDFs
            PDFs(i,:) = exp(fnval(spmak(knots,pchain(i,:)),x));
        else
            PDFs(i,:) = fnval(spmak(knots,pchain(i,:)),x);
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
    
    %initialize histogram grid
    histgrid = zeros(y_res,x_res);
    
    %determine edges of the bin spacing in the y (log P) direction
    minP = min(min(PDFs));
    maxP = max(max(PDFs));
    y_edges = [minP:(maxP-minP)/y_res:maxP];
    y = [y_edges(2:end)-(y_edges(2)-y_edges(1))/2];
    yzero = find(y==min(abs(y)));
    
    %calculate the number of PDFs passing through each grid cell
    h = parfor_progressbar(length(x),'Calculating...');
    parfor i = 1:length(x)
        histgrid(:,i) = histcounts(PDFs(:,i),y_edges)';
        h.iterate();
    end
    close(h);
    %limit histogram bin values so that a few cells with extremely high
    %values don't affect the color plotting for the whole figure.  Also
    %ensure that function values near 0 do not dominate the color scale if
    %the p scale is linear.
    if ~logPDFs
        histgridnonzero = histgrid;
        histgridnonzero(:,yzero)=[];
        histgridvals = reshape(histgridnonzero,[],1);
    else
        histgridvals = histgrid(:);
    end
    hist95 = quantile(histgridvals(find(histgridvals~=0&histgridvals~=max(histgridvals))),.95);
    histgrid(find(histgrid>hist95)) = hist95;
    
    [X,Y] = meshgrid(x,y);
    
    figure(fignum);
    clf;
    %alters colormap so that bins that zero PDFs pass through plot as white
    cmap = colormap;
    cmap(1,:) = [1 1 1];
    colormap(cmap);
    
    %output histogram values showing concentration of PDF curves of the
    %PME.
    surf(exp(X),Y,histgrid,'EdgeColor','none');
    set(gca, 'XScale', 'log');
    ax = gca;
    ax.XLim = [exp(x_min) exp(x_max)];
    yax = ax.YLim;
    xlabel('Age (Ma)');
    if logPDFs
        ylabel('ln p');
    else
        ylabel('p');
    end
    view(2);
    
    %plot the maximum likelihood model PDF
    if(showMLM==1)
        hold on;
        %max. likelihood model can be output onto a new figure
        if newfig
            figure(fignum+1);
            clf
            ax = gca;
            if logPDFs
                h=plot3(exp(x),fnval(spmak(knots,pchain(1,:)),x)+log((x_max-x_min)/(exp(x_max)-exp(x_min))),hist95+1*ones(1,length(x)),'r-','LineWidth',1);
            else
                h=plot3(exp(x),exp(fnval(spmak(knots,pchain(1,:)),x))*(x_max-x_min)/(exp(x_max)-exp(x_min)),hist95+1*ones(1,length(x)),'r-','LineWidth',1);
            end
            set(gca, 'XScale', 'log');
            ax.YLim = yax;
            ax.XLim = [exp(x_min) exp(x_max)];
            view(2)
            xlabel('Age (Ma)');
            if logPDFs
                ylabel('ln p');
            else
                ylabel('p');
            end
        %or can be plotted onto the same figure as the PME.
        else
            if logPDFs
                h=plot3(exp(x),fnval(spmak(knots,pchain(1,:)),x)+log((x_max-x_min)/(exp(x_max)-exp(x_min))),hist95+1*ones(1,length(x)),'r-','LineWidth',1);
            else
                h=plot3(exp(x),exp(fnval(spmak(knots,pchain(1,:)),x))*(x_max-x_min)/(exp(x_max)-exp(x_min)),hist95+1*ones(1,length(x)),'r-','LineWidth',1);
            end
        end
    end
    %plot dots corresponding to measured zircon ages
    if showdot        
        hold on;
        ages = csvread(strcat(folderpath,name,'.csv'));
        %eliminate index column in age data, if it exists.
        if size(ages,2)==3
            ages = ages(:,[2 3]);
        end
        %make marker size contingent on number of datapoints being shown
        if size(ages,1)>300
            msize = 5;
        else
            msize = 10;
        end
        if newfig
            figure(fignum+2);
            clf
            plot(ages(:,1),rand(1,size(ages,1)),'k.','MarkerSize',msize);
            set(gca, 'XScale', 'log');
            axis([exp(x_min) exp(x_max) 0 1]);
            xlabel('Age (Ma)');
        else
            plot(ages(:,1),minP-((maxP-minP)/20)*rand(1,size(ages,1)),'k.','MarkerSize',msize);
        end
    end
    %display title
    if showtitle
        figure(fignum);
        title(plottitle,'Interpreter', 'none');
    end
    z=0;
    
end