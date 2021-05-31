function plotCM =  plotColor(nPlot,varargin)

    switch nargin
        case 1
            cm = colormap;
        case 2
            cm = colormap(varargin{1});
    end

    ncm = size(cm,1);
    indPlot = 1:nPlot;
    plotCM = cm(round((indPlot-1)/(nPlot-1)*(ncm-1)+1),:);
    
end
















