function hasError = dispNaNInf(mat)

    varName = inputname(1);
    hasError = 0;
    if any(isnan(mat(:)))
        charMat = 'NaN';
        indNaN = find(isnan(mat(:)));
        nNaN = length(indNaN);
        warning([num2str(nNaN) ' elements are ' charMat ' in ' varName]);
        hasError = 1;
    end
    if any(isinf(mat(:)))
        charMat = 'Inf';
        indInf = find(isinf(mat(:)));
        nInf = length(indInf);
        warning([num2str(nInf) ' elements are ' charMat ' in ' varName]);
        hasError = 1;
    end    
    

end