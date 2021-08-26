
%% set up the points
clear all; close all; clc;

nPoints = 1000;
vectorLength = 128;
noiseScale = 0.2;
points1 = 2*rand(nPoints, vectorLength) - 1;
% Noise generated from a uniform random distribution.
% noise = (2*rand(nPoints, vectorLength) - 1) * noiseScale;
% Noise generated from a normal distribution
noise = randn(nPoints, vectorLength) * noiseScale;
points2 = points1 + noise;

% Normalize the vectors so their absolute values are all 1 (l2 norms are 1)
for i = 1:nPoints
    points1(i,:) = points1(i,:) / norm(points1(i,:), 2);
    points2(i,:) = points2(i,:) / norm(points2(i,:), 2);
end

%% matching using naive comparisons
naiveMatchingStart = tic;

matching_naive = zeros(nPoints, 1);

for i = 1:nPoints
    bestMatchDist = 1e10;
    match = -1;
    for j = 1:nPoints
        % Take euclidean distance between the vector
        diff = sum((points1(i,:) - points2(j,:)) .^ 2);
        % If euclidean distance is better than the best match,
        % we have a new best match
        if (diff < bestMatchDist)
            bestMatchDist = diff;
            match = j;
        end
    end
    matching_naive(i) = match;
end
naiveMatchingTime = toc(naiveMatchingStart);

%% Contstruct LSH table (one LSH table)
lshStart = tic;

nTables = 8;

% Initialize random hyperplanes
nPlanesPerTable = round(log2(nPoints));
hyperplanes = 2*rand(nTables*nPlanesPerTable, vectorLength) - 1;
maxHPLength = 1;

% Initialize LSH tables
indexGroupMap = zeros(nTables, nPoints);
nBoxes = 2^nPlanesPerTable;
groupSizeMap = zeros(nTables, nBoxes);
groupIndexMap = zeros(nTables, nBoxes);
groupIndexMapTails = zeros(nTables, nBoxes);

% array of points grouped by what box they got hashed into
groupArray = zeros(nPoints, 1);

% - Put points1 into our lsh table -

% Calculate hash values, keep track of sizes of each group
for table = 1:nTables
    for i = 1:nPoints
        point = points1(i,:);
        hashcode = 0;
        for j = 1:nPlanesPerTable
            hplane = hyperplanes((table-1)*nPlanesPerTable + j,:)';
            if (point*hplane > 0)
                hashcode = bitor(hashcode, bitshift(1, j-1));
            end
        end
        hashcode = hashcode+1; % Because matlab is weird with indexing
        indexGroupMap(table, i) = hashcode;
        groupSizeMap(table, hashcode) = groupSizeMap(table, hashcode)+1;
    end

    % Prepare the index map
    cnt = 1;
    for i = 1:nBoxes
        groupIndexMap(table, i) = cnt;
        groupIndexMapTails(table, i) = groupIndexMap(table, i);
        cnt = cnt + groupSizeMap(table, i);
    end


    for i = 1:nPoints
        % What group did the point get hashed into
        hashcode = indexGroupMap(table, i);
        % What is the start index of that group in groupArray
        idx = groupIndexMapTails(table, hashcode);

        % Increment the index tail mapping
        groupIndexMapTails(table, hashcode) = idx + 1;

        groupArray(table, idx) = i;
    end
end

%% Match points2 with points1 using LSH hash table

% Calculating hash codes happens separately from searching,
% this makes the program easier to parallelize at some later point

bestMatchDists = zeros(nPoints) + 1e10;
matching_lsh = zeros(nPoints, 1);

% Calculate hash values
for table = 1:nTables
    for i = 1:nPoints
        point = points2(i,:);
        hashcode = 0;
        for j = 1:nPlanesPerTable
            hplane = hyperplanes((table-1)*nPlanesPerTable + j,:)';
            if (point*hplane > 0)
                hashcode = bitor(hashcode, bitshift(1, j-1));
            end
        end
        hashcode = hashcode+1; % Because matlab is weird with indexing
        indexGroupMap(table, i) = hashcode;
    end

    for i = 1:nPoints
        % Find table and elements in points1 to match with
        changed = 0;
        bestMatchDist = bestMatchDists(i);
        hashcode = indexGroupMap(table, i);
        size = groupSizeMap(table, hashcode);
        startIdx = groupIndexMap(table, hashcode);

        % Match the points
        for j = startIdx:(startIdx+size-1)
            % Take euclidean distance between the vector
            idx = groupArray(table, j);
            diff = sum((points2(i,:) - points1(idx,:)) .^ 2);
            % If euclidean distance is better than the best match,
            % we have a new best match
            if (diff < bestMatchDist)
                bestMatchDist = diff;
                match = idx;
                changed = 1;
            end
        end
        if (changed == 1)
%             fprintf("prev value: %d    new value: %d\n", matching_lsh(i), match)
            matching_lsh(i) =  match;
            bestMatchDists(i) = bestMatchDist;
        end
    end
end

lshTime = toc(lshStart);


%% Deal with time

wrong = 0;
for i = 1:nPoints
    if (matching_lsh(i) ~= matching_naive(i))
        wrong = wrong + 1;
    end
end
correct = nPoints - wrong;
correctRatio = correct/nPoints;

fprintf("N = %d\n", nPoints);
fprintf("Matching time using naive method:      %f\n", naiveMatchingTime);
fprintf("Total matching time using LSH:         %f\n", lshTime);
fprintf("\n\n");
fprintf("Time ratio (lsh time / naive time): %f\n", (lshTime/naiveMatchingTime));
fprintf("LSH correct ratio:                  %f\n", correctRatio);


