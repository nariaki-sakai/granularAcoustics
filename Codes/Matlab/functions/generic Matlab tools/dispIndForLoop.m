function dispIndForLoop(ind,indMax,varargin)
    if nargin == 2
        step = 10^(ceil(log10(indMax/10)));
    else
        step = varargin{1};
    end
    if mod(ind,step) == 0
        disp([num2str(ind) ' / ' num2str(indMax)]);
    end
end