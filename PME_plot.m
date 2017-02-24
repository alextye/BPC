function z = PME_plot(folderpath,name,x_min,x_max,x_res,y_res,fignum)
    
%OUTPUT
%y is a dummy variable

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
    for i = 1:size(PDFs,1)
        PDFs(i,:) = fnval(spmak(knots,pchain(i,:)),x);
    end
    
    %initialize histogram grid
    histgrid = zeros(y_res,x_res);
    
    %determine edges of the bin spacing in the y (log P) direction
    minP = min(min(PDFs));
    maxP = max(max(PDFs));
    y_edges = [minP:(maxP-minP)/y_res:maxP];
    y = [y_edges(2:end)-(y_edges(2)-y_edges(1))/2];
    
    for i = 1:length(x)
        %keyboard
        histgrid(:,i) = histcounts(PDFs(:,i),y_edges)';
    end
    
    %alters colormap so that no data plots as white
    cmap = colormap;
    cmap(1,:) = [1 1 1];
    colormap(cmap);
    
    [X,Y] = meshgrid(x,y);
    
    figure(fignum);
    
    mesh(X,Y,histgrid,'EdgeColor','none');
    view(2);
    keyboard
    z=0;
    %imagesc([x(1) x(end)],[y(end) y(1)],histgrid);
end