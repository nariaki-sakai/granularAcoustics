function rank = detectorPos2ind(detectorPos,Ltot)

    %%
    
    nCell = size(detectorPos);
    lenPos = zeros(nCell);
    for ind = 1:nCell
        lenPos(ind) = length(detectorPos{ind});
    end
    
    subPos = cell2mat(detectorPos);
    indPos = sub2ind(Ltot,subPos(:,1),subPos(:,2));
    [indPos,~,ic] = unique(indPos);
    [~, rank] = sort(indPos,'ascend');
    rank = rank(ic);
    rank = mat2cell(rank,lenPos);
    
    
end


















