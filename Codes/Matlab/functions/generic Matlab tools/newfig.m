function h = newfig(nFig,varargin)


    figure(nFig);
    close(nFig);
    h = figure(nFig);
%     set(gca,'plotBoxAspectRatio',[1 1 1]);
    hold all;
    
    if nargin == 2
        windowSize = varargin{1};
        formatfig(nFig,windowSize);
    else
        formatfig(nFig);
    end
    

end








