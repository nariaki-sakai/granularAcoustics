function [f, yfft] = fourierTransform(t,y,nt,optNorm)

    %%
    dt = t(2) - t(1);
%     nt = length(t);

    if isrow(y)
        y = y';
    end
    yfft = fft(y,nt,1);

    fs = 1/(dt*nt);
    if mod(nt,2)
        f = fs*(0:(nt-1)/2);
        yfft = yfft(1:(nt-1)/2+1,:);
        yfft(2:end,:) = 2*yfft(2:end,:);
    else
        f = fs*(0:nt/2);
        yfft = yfft(1:nt/2+1,:);
        yfft(2:end,:) = 2*yfft(2:end,:);
    end
    
    switch optNorm
        case 'power'
            yfft = yfft/nt;
        case 'energy'
            yfft = yfft*dt;
        otherwise
            error('optNorm has wrong attribute');
    end
%     yfft(2:nt/2+1) = 2*yfft(2:nt/2+1);

end












































