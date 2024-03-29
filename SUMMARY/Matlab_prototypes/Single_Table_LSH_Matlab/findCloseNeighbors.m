function [N, hashcodes, combDists] = ...
    findCloseNeighbors(sqrdDists, threshold, nPlanes, hashcode)

    [setLen, combMasks, combDists] = findNeighborMasks(sqrdDists, threshold, nPlanes);
    N = setLen+1;
    
    hashcodes = zeros(N);
    hashcodes(1) = hashcode;
    for i = 1:setLen
        hashcodes(i+1) = bitxor(hashcode, combMasks{i});
    end
end
