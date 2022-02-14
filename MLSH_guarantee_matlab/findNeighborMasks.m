function [setLen, combMasks, combDists] = ...
    findNeighborMasks(sqrdDists, threshold, nPlanes)

    % threshold = 6;
    % nPlanes = 5;
    % sqrdDists = [3, 1, 4, 8, 2];

    % Masks and dists for combination of planes
    combMasks = cell(1,20);
    combDists = cell(1,20);

    % Add single values below the threshold to lists
    setLen = 0;
    for i = 1:nPlanes
        if (sqrdDists(i) < threshold)
            setLen = setLen + 1;
            combMasks{setLen} = bitshift(1,i-1);
            combDists{setLen} = sqrdDists(i);
        end
    end

    % For every value starting with the last, try to combine the value
    % with all other values in the list
    underThreshLen = setLen;
    for k = (underThreshLen-1):-1:1
        tmp = setLen;
        for i = (k+1):setLen
            if (combDists{k} < threshold)
                dist = combDists{k} + combDists{i};
                if (dist < threshold)
                    tmp = tmp + 1;
                    combMasks{tmp} = bitor(combMasks{k}, combMasks{i});
                    combDists{tmp} = dist;
                end
            end
        end
        setLen = tmp;
    end

    % for i = 1:setLen
    %     fprintf("hashcode = %s,  distance = %d\n", ...
    %         dec2bin(combMasks{i}, nPlanes), combDists{i});
    % end

end