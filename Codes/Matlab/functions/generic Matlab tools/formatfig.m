function formatfig(varargin)

    if nargin == 1
        if length(varargin{1}) == 1
            nf = varargin{1};
            windowSize = [0 0 300 300];
        else
            nf = gcf;
            windowSize = varargin{1};
            windowSize = [0 0 windowSize];
        end
    else
        nf = varargin{1};
        windowSize = varargin{2};
        windowSize = [0 0 windowSize];
    end
    figure(nf);
    screenSize = get(0,'ScreenSize');
    
    posWin = get(gcf,'Position');
%     windowSize = [randi(screenSize(3)-windowSize(3)) randi(screenSize(4)-windowSize(4)) windowSize(3:4)];
    windowSize = [posWin(1:2) windowSize(3:4)];
    set(gcf,'Position',windowSize);

    set(gcf,'defaultLegendAutoUpdate','off');
    
    hg = findall(gcf,'type','axes');
    if length(hg) < 2
        for indg = 1:length(hg)
            %pbaspect(hg(indg),[1 1 1]);
            
            set(hg(indg),'fontsize',16);
            box(hg(indg),'on');
        end
    else
        for indg = 1:length(hg)
            set(hg(indg),'fontsize',12);
            box(hg(indg),'on');
        end
    end
end