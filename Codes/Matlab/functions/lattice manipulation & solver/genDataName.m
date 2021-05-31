function charData = genDataName(L,p,dk,strain,a,indSample,f,forcing)

    %charData = genDataName(L,p,dk,strain,a,indSample,f,forcing)
    %   generates a string of character from the parameters of
    %   the numerical simulation. Useful for e.g. producing a formatted
    %   file name for data or figures
    %
    %   INPUT: see loadLattice, inputs are almost the same
    %
    %   OUTPUT:
    %       charData: formatted string of character. See the code

    if f > 0
        T = 1/f;
        if T >= 1
            charData = sprintf('L%dp%d%%dk%d%%strain%d%%springL%d%%sample%dT%d%s',L(1),100*p,100*dk,100*strain,100*a,indSample,T,forcing);
        else
            charData = sprintf('L%dp%d%%dk%d%%strain%d%%springL%d%%sample%df%d%s',L(1),100*p,100*dk,100*strain,100*a,indSample,f,forcing);
        end
    else
        charData = sprintf('L%dp%d%%dk%d%%strain%d%%springL%d%%sample%d%s',L(1),100*p,100*dk,100*strain,100*a,indSample,forcing);
    end
    charData = strrep(charData,'.','_');
    
end