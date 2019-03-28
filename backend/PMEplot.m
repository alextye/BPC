function z = PMEplot(folderpath,name,x_min,x_max,x_res,y_res,fignum,showMLM,showdot,logPDFs,thin,linage,hybrid,varargin)
    
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

    %read knot locations for spline
    %N_coefs=50; %by default
    %knots = [x_min x_min x_min:(x_max-x_min)/(N_coefs-2):x_max x_max x_max];
    data = csvread(strcat(folderpath,'chains/',name,'chain.csv'),0);
    knots = data(1,:);
    
    %if min/max limits are outside the modeled range, they will be set to
    %lie at the limits of the modeled range.
    x_min = max(x_min,knots(1));
    x_max = min(x_max,knots(end));
    %read pchain
    pchain = data([2:end],[1:end-3]);
    %keyboard
    
    %thin (select a subset of models) for plotting if user has specified to
    %do so.
    if thin > 0 & thin < size(pchain,1)
        %preserve MLM at head of pchain
        pchain1 = pchain(1,:);
        
        idx = randperm(size(pchain,1)-1,thin);
        pchain = pchain(idx,:);
        pchain(1,:) = pchain1;
    end
    
    if hybrid
        x_break = min(x_max,log(1000));
        x1 = [x_min+(x_break-x_min)/(2*x_res):(x_break-x_min)/x_res:x_break-(x_break-x_min)/(2*x_res)];
        if x_max > x_break
            x2 = log([exp(x_break)+(exp(x_max)-exp(x_break))/(2*x_res):(exp(x_max)-exp(x_break))/x_res:exp(x_max)-(exp(x_max)-exp(x_break))/(2*x_res)]);
        else
            x2 = [];
        end
        x = [x1 x2];
        PDFs = zeros(size(pchain,1),length(x));
    elseif linage
        x = log([exp(x_min)+(exp(x_max)-exp(x_min))/(2*x_res):(exp(x_max)-exp(x_min))/x_res:exp(x_max)-(exp(x_max)-exp(x_min))/(2*x_res)]);
        PDFs = zeros(size(pchain,1),x_res);
    else
        %set x values at which to evaluate each PDF, dictated by x_min,
        %x_max, and x_res
        x = [x_min+(x_max-x_min)/(2*x_res):(x_max-x_min)/x_res:x_max-(x_max-x_min)/(2*x_res)];
        PDFs = zeros(size(pchain,1),x_res);
    end
    
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
        if hybrid
            PDFs(:,find(x>x_break)) = PDFs(:,find(x>x_break))./(exp(x(find(x>x_break)))/1000);
        elseif linage
            PDFs = PDFs./exp(x);
        else
            PDFs = PDFs * (knots(end)-knots(1))/(exp(knots(end))-exp(knots(1)));
        end
    else
        if hybrid
            PDFs(:,find(x>x_break)) = PDFs(:,find(x>x_break))-log(exp(x(find(x>x_break)))/1000);
        elseif linage
            PDFs = PDFs-x;
        else
            PDFs = PDFs + log((knots(end)-knots(1))/(exp(knots(end))-exp(knots(1))));
        end
    end
    
%    keyboard
    
    %initialize histogram grid
    histgrid = zeros(y_res,length(x));
    
    %determine edges of the bin spacing in the y (log P) direction
    %minP = min(min(PDFs));
    minP = quantile(min(PDFs),.01);
    %maxP = max(max(PDFs));
    maxP = quantile(max(PDFs),.99);
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
    
    %STILL NEED TO FIGURE OUT HOW TO PLOT 2 PLOTS FOR THE HYBRID THING
    figure(fignum);
    clf;
    %alters colormap so that bins that zero PDFs pass through plot as white
    cmap = colormap;
    cmap(1,:) = [1 1 1];
    colormap(cmap);

    %output histogram values showing concentration of PDF curves of the
    %PME.
    if ~hybrid
        surf(exp(X),Y,histgrid,'EdgeColor','none');
    else
        surf(exp(X(:,[1:x_res])),Y(:,[1:x_res]),histgrid(:,[1:x_res]),'EdgeColor','none');
    end
    if ~linage
        set(gca, 'XScale', 'log');
    end
    ax = gca;
    
    if ~hybrid
        ax.XLim = [exp(x_min) exp(x_max)];
    else
        ax.XLim = [exp(x_min) exp(x_break)]; 
    end
    
    yax = ax.YLim;
    xlabel('Age (Ma)');
    if logPDFs
        ylabel('ln p');
    else
        ylabel('p');
    end
    view(2);
    grid off
    %keyboard
    if hybrid
        figure(fignum+1);
        clf;
        %alters colormap so that bins that zero PDFs pass through plot as white
        cmap = colormap;
        cmap(1,:) = [1 1 1];
        colormap(cmap);

        %output histogram values showing concentration of PDF curves of the
        %PME.
        surf(exp(X(:,[x_res+1:end])),Y(:,[x_res+1:end]),histgrid(:,[x_res+1:end]),'EdgeColor','none');
        %set(gca, 'XScale', 'log');
        ax = gca;
        ax.XLim = [exp(x_break) exp(x_max)];
        yax = ax.YLim;
        xlabel('Age (Ma)');
        if logPDFs
            ylabel('ln p');
        else
            ylabel('p');
        end
        view(2);
        grid off
    end
    
    %plot the maximum likelihood model PDF
    if(showMLM==1)
        hold on;
        %max. likelihood model can be output onto a new figure
        if newfig
            if ~hybrid
                figure(fignum+1);
            else
                figure(fignum+2);
            end
            clf
            ax = gca;
            h=plot3(exp(x),PDFs(1,:),hist95+1*ones(1,length(x)),'r-','LineWidth',1);
            %keyboard
            if ~linage
                set(gca, 'XScale', 'log');
            end
            ax.YLim = yax;
            if ~hybrid
                ax.XLim = [exp(x_min) exp(x_max)];
            else
                ax.XLim = [exp(x_min) exp(x_break)];
            end
            view(2)
            grid off
            xlabel('Age (Ma)');
            if logPDFs
                ylabel('ln p');
            else
                ylabel('p');
            end
            if hybrid
                figure(fignum+3);
                clf
                ax = gca;
                h=plot3(exp(x),PDFs(1,:),hist95+1*ones(1,length(x)),'r-','LineWidth',1);
                %keyboard
                ax.YLim = yax;
                ax.XLim = [exp(x_break) exp(x_max)];
                view(2)
                grid off
                xlabel('Age (Ma)');
                if logPDFs
                    ylabel('ln p');
                else
                    ylabel('p');
                end
            end
        %or can be plotted onto the same figure as the PME.
        else
            if ~hybrid
                %keyboard
                h=plot3(exp(x),PDFs(1,:),hist95+1*ones(1,length(x)),'r-','LineWidth',1);
            else
                h=plot3(exp(x),PDFs(1,[1:x_res]),hist95+1*ones(1,length(x)),'r-','LineWidth',1);
                figure(fignum+1)
                h=plot3(exp(x),PDFs(1,[x_res+1:end]),hist95+1*ones(1,length(x)),'r-','LineWidth',1);
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
            if ~hybrid
                figure(fignum+2);
            else
                figure(fignum+4);
            end
            clf
            plot(ages(:,1),rand(1,size(ages,1)),'k.','MarkerSize',msize);
            if ~linage
                set(gca, 'XScale', 'log');
            end
            axis([exp(x_min) exp(x_max) 0 1]);
            xlabel('Age (Ma)');
            if hybrid
                axis([exp(x_min) exp(x_break) 0 1]);
                figure(fignum+5);
                clf
                plot(ages(:,1),rand(1,size(ages,1)),'k.','MarkerSize',msize);
                axis([exp(x_break) exp(x_max) 0 1]);
                xlabel('Age (Ma)');
            end
        else
            if ~hybrid
                plot(ages(:,1),minP-((maxP-minP)/20)*rand(1,size(ages,1)),'k.','MarkerSize',msize);
            else
                figure(fignum)
                ages_lessbreak = ages(find(ages(:,1)<exp(x_break)),:);
                ages_gtrbreak = ages(find(ages(:,1)>exp(x_break)),:);
                plot(ages_lessbreak(:,1),minP-((maxP-minP)/20)*rand(1,size(ages_lessbreak,1)),'k.','MarkerSize',msize);
                figure(fignum+1)
                plot(ages_gtrbreak(:,1),minP-((maxP-minP)/20)*rand(1,size(ages_gtrbreak,1)),'k.','MarkerSize',msize);
            end
        end
    end
    %display title
    if showtitle
        figure(fignum);
        title(plottitle,'Interpreter', 'none');
    end
    z=0;
    
end