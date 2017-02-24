function z = PMEplot(folderpath,name,x_min,x_max,x_res,y_res,fignum,varargin)
    
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
%y-coordinates.

%fignum--this is the number of the figure to which to output the display

%varargin--a value of zero or 1 may be optionally input to suppress the
%plotting of the maximum likelihood model line in the plot.  Inputting 1
%will cause the line not to be plotted and inputting 0 will cause it to be
%plotted.

    if size(varargin,2)>0
        supline = varargin{1};
    else
        supline = 0;
    end

    %log transform x_min and x_max
    x_min = log(x_min);
    x_max = log(x_max);

    %set knot locations for spline
    N_coefs=50; %by default
    knots = [x_min x_min x_min:(x_max-x_min)/(N_coefs-2):x_max x_max x_max];
    %read pchain
    pchain = csvread(strcat(folderpath,'chains/',name,'chain.csv'),0);
    %set x values at which to evaluate each PDF, dictated by x_min,
    %x_max, and x_res
    x = [x_min+(x_max-x_min)/(2*x_res):(x_max-x_min)/x_res:x_max-(x_max-x_min)/(2*x_res)];
    PDFs = zeros(size(pchain,1),x_res);
    
    %generate logPDF function values at indicated x-values for all
    %probability chains
    parfor i = 1:size(PDFs,1)
        PDFs(i,:) = fnval(spmak(knots,pchain(i,:)),x);
    end
    
    %initialize histogram grid
    histgrid = zeros(y_res,x_res);
    
    %determine edges of the bin spacing in the y (log P) direction
    minP = min(min(PDFs));
    maxP = max(max(PDFs));
    y_edges = [minP:(maxP-minP)/y_res:maxP];
    y = [y_edges(2:end)-(y_edges(2)-y_edges(1))/2];
    
    %calculate the number of PDFs passing through each grid cell
    parfor i = 1:length(x)
        histgrid(:,i) = histcounts(PDFs(:,i),y_edges)';
    end
    %limit histogram bin values so that a few cells with extremely high
    %values don't affect the color plotting for the whole figure
    histgridvals = reshape(histgrid,[],1);
    hist99 = quantile(histgridvals(find(histgridvals~=0&histgridvals~=max(histgridvals))),.99);
    histgrid(find(histgrid>hist99)) = hist99;
    
    [X,Y] = meshgrid(x,y);
    
    figure(fignum);
    %alters colormap so that bins that zero PDFs pass through plot as white
    cmap = colormap;
    cmap(1,:) = [1 1 1];
    colormap(cmap);
    
    %output histogram values with log scales.
    surf(exp(X),Y,histgrid,'EdgeColor','none');
    set(gca, 'XScale', 'log');
    ax = gca;
    ax.XLim = [exp(x_min) exp(x_max)];
    xlabel('Age (Ma)');
    ylabel('ln p');
    view(2);
    %plot the maximum likelihood model PDF
    if(supline==0)
        hold on;
        h=plot3(exp(x),fnval(spmak(knots,pchain(1,:)),x),hist99+1*ones(1,length(x)),'r-','LineWidth',1);        
    end
        
    z=0;
end