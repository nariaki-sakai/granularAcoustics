function modeRange = guideModeIdxRange(f,L,strain)
    
    %%
    omega0 = 2*pi*sqrt(strain/(1+strain));
    L = L(1);
    
    omega = 2*pi*f;
    if omega < 2*omega0
        nmin = 1;
        nmax = 2*(L-1)/pi*asin(omega/(2*omega0));
        if mod(nmax,1) == 0
            nmax = nmax - 1;
        end
    elseif omega < 2*sqrt(2)*omega0
        nmin = 2*(L-1)/pi*asin(sqrt((omega/(2*omega0))^2 - 1));
        nmax = L-2;
    else
        error('Excitation frequency out of max limit');
    end
    
    modeRange = [nmin,nmax];

end
































