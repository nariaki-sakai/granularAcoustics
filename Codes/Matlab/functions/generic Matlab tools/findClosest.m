function [ind, value] = findClosest(xTarget,xVector)
    if length(xTarget) == 1
        [~, ind]  = min(abs(xVector - xTarget));
        value = xVector(ind);
    else
        disp('ERROR: second input should be a single variable');
    end
end