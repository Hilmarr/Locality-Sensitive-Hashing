%% set up the points
clear all; close all; clc;

% maxHPLengthArr = [0, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3];
maxHPLengthArr = linspace(0, 0.003, 10);

nPoints = 5000;
vectorLength = 128;
points1 = 2*rand(nPoints, vectorLength) - 1;

% Normalize the vectors so their absolute values are all 1 (l2 norms are 1)
for i = 1:nPoints
    points1(i,:) = points1(i,:) / norm(points1(i,:), 2);
end


%% Contstruct LSH table (one LSH table). Test number of duplicates
nDuplicatesArr = zeros(length(maxHPLengthArr));
largestBoxArr = zeros(length(maxHPLengthArr));

for curIter = 1:length(maxHPLengthArr)
    
    maxHPLength = maxHPLengthArr(curIter);

    % Initialize random hyperplanes
    nPlanes = round(log2(nPoints));
    hyperplanes = 2*rand(nPlanes, vectorLength) - 1;
    hyperplaneLengths = maxHPLength * rand(1, nPlanes);
    
    % Give the hyperplanes the length of somewhere between 0 and maxHPLength
    % Not needed if all hyperplanes go through the origin
    for i = 1:nPlanes
        hyperplanes(i,:) = hyperplanes(i,:) / norm(hyperplanes(i,:), 2);
    end

    % Initialize LSH tables
    indexGroupMap = zeros(nPoints, 1);
    nBoxes = 2^nPlanes;
    groupSizeMap = zeros(nBoxes, 1);

    % array of points grouped by what box they got hashed into
    groupArray = zeros(nPoints, 1);

    % - Put points1 into our lsh table -

    % Calculate hash values, keep track of sizes of each group
    for i = 1:nPoints
        point = points1(i,:);
        hashcode = 0;
        for j = 1:nPlanes
            hplane = hyperplanes(j,:)';
            if (point*hplane > hyperplaneLengths)
                hashcode = bitor(hashcode, bitshift(1, j-1));
            end
        end
        hashcode = hashcode+1; % Because matlab is weird with indexing
        groupSizeMap(hashcode) = groupSizeMap(hashcode)+1;
    end

    % Count duplicates in indexGroupMap
    % i.e. the number of points that gets mapped to a box where there already
    %      exists some points
    % (in order to check whether using hyperplanes that go through the origin
    %  is better than more random hyperplanes, should probably plot number
    %  of duplicates given the length of the hyperplanes)
    nDuplicatesArr(curIter) = sum(groupSizeMap(groupSizeMap > 1)-1);
    % See the maximum amount of points that get mapped to the same box
    largestBoxArr(curIter) = max(groupSizeMap);

end

figure
plot(maxHPLengthArr, nDuplicatesArr);
title(sprintf("Number of duplicates as a function of hyperplane length, N = %d", nPoints))
xlabel("max hyperplane length");
ylabel("number of duplicates");
% set(gca, 'XScale', 'log')

figure
plot(maxHPLengthArr, largestBoxArr);
title(sprintf("Number of duplicates as a function of hyperplane length, N = %d", nPoints))
xlabel("max hyperplane length");
ylabel("number of duplicates");
