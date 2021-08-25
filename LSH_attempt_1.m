
%% set up the points
clear all; close all; clc;

nPoints = 1000;
vectorLength = 128;
noiseScale = 0.01;
points1 = 2*rand(nPoints, vectorLength) - 1;
% Noise generated with uniform random distribution.
% Might later try it with a normal distribution, but this is good for now.
noise = (2*rand(nPoints, vectorLength) - 1) * noiseScale;
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
    bestMatch = 1e10;
    match = -1;
    for j = 1:nPoints
        % Take euclidean distance between the vector
        diff = sum((points1(i,:) - points2(j,:)) .^ 2);
        % If euclidean distance is better than the best match,
        % we have a new best match
        if (diff < bestMatch)
            bestMatch = diff;
            match = j;
        end
    end
    matching_naive(i) = match;
end
naiveMatchingTime = toc(naiveMatchingStart);

%% Contstruct LSH table (one LSH table)
lshStart = tic;

% Initialize random hyperplanes
nPlanes = round(log2(nPoints));
hyperplanes = 2*rand(nPlanes, vectorLength) - 1;
maxHPLength = 1;

% % Give the hyperplanes the length of somewhere between 0 and maxHPLength
% % Not needed if asll hyperplanes go through the origin
% for i = 1:nPlanes
%     hyperplanes(i,:) = hyperplanes(i,:) ...
%         * (maxHPLength * rand(1) / norm(hyperplanes(i,:), 2));
% end

% Initialize LSH tables
indexGroupMap = zeros(nPoints, 1);
nBoxes = 2^nPlanes;
groupSizeMap = zeros(nBoxes, 1);
groupIndexMap = zeros(nBoxes, 1);
groupIndexMapTails = zeros(nBoxes, 1);

% array of points grouped by what box they got hashed into
groupArray = zeros(nPoints, 1);

% - Put points1 into our lsh table -

% Calculate hash values, keep track of sizes of each group
hashValuesStart1 = tic;
for i = 1:nPoints
    point = points1(i,:);
    hashcode = 0;
    for j = 1:nPlanes
        hplane = hyperplanes(j,:)';
%         if (point*hplane > hplane'*hplane)
%             hashcode = bitor(hashcode, bitshift(1, j-1));
%         end
        if (point*hplane > 0)
            hashcode = bitor(hashcode, bitshift(1, j-1));
        end
    end
    hashcode = hashcode+1; % Because matlab is weird with indexing
    indexGroupMap(i) = hashcode;
    groupSizeMap(hashcode) = groupSizeMap(hashcode)+1;
end
hashValuesTime1 = toc(hashValuesStart1);

% % LATER: count duplicates in indexGroupMap
% % i.e. the number of points that gets mapped to a box where there already
% %      exists some points
% % (in order to check whether using hyperplanes that go through the origin
% %  is better than more random hyperplanes, should probably plot number
% %  of duplicates given the length of the hyperplanes)
% nDuplicates = sum(groupSizeMap(groupSizeMap > 1)-1);
% % See the maximum amount of points that get mapped to the same box
% largestBox = max(groupSizeMap);

% Prepare the index map
cnt = 1;
for i = 1:nBoxes
    groupIndexMap(i) = cnt;
    groupIndexMapTails(i) = groupIndexMap(i);
    cnt = cnt + groupSizeMap(i);
end

for i = 1:nPoints
    % What group did the point get hashed into
    hashcode = indexGroupMap(i);
    % What is the start index of that group in groupArray
    idx = groupIndexMapTails(hashcode);
    
%     % Add to table
%     if (idx > nPoints || idx < 1)
%         idx
%     end
%     if (groupArray(idx) > 0)
%         fprintf("HEY!!! groupArray(idx) == %d, now putting in %d\n", groupArray(idx), i);
%     end
    
    % Increment the index tail mapping
    groupIndexMapTails(hashcode) = idx + 1;
    
    groupArray(idx) = i;
    
end

tablesCreatedTime = toc(lshStart);

%% Match points2 with points1 using LSH hash table
matchingStart = tic;

% Calculating hash codes happens separately from searching,
% this makes the program easier to parallelize at some later point

% Calculate hash values
hashValuesStart2 = tic;
for i = 1:nPoints
    point = points2(i,:);
    hashcode = 0;
    for j = 1:nPlanes
        hplane = hyperplanes(j,:)';
        if (point*hplane > 0)
            hashcode = bitor(hashcode, bitshift(1, j-1));
        end
    end
    hashcode = hashcode+1; % Because matlab is weird with indexing
    indexGroupMap(i) = hashcode;
end
hashValuesTime2 = toc(hashValuesStart2);

% Match the points using lsh
matching_lsh = zeros(nPoints, 1);

for i = 1:nPoints
    % Find table and elements in points1 to match with
    bestMatch = 1e10;
    hashcode = indexGroupMap(i);
    size = groupSizeMap(hashcode);
    startIdx = groupIndexMap(hashcode);
    
    % Match the points
    for j = startIdx:(startIdx+size-1)
        % Take euclidean distance between the vector
        idx = groupArray(j);
        diff = sum((points2(i,:) - points1(idx,:)) .^ 2);
        % If euclidean distance is better than the best match,
        % we have a new best match
        if (diff < bestMatch)
            bestMatch = diff;
            match = idx;
        end
    end
    if (bestMatch ~= 1e10)
        matching_lsh(i) =  match;
    else
        matching_lsh(i) = -1;
    end
end

matchingTime = toc(matchingStart);
lshTime = toc(lshStart);
hashValuesTime = hashValuesTime1 + hashValuesTime2;


% Deal with time

naiveMatchingTime
tablesCreatedTime
matchingTime
lshTime
hashValuesTime1
hashValuesTime2
hashValuesTime


wrong = 0;
for i = 1:nPoints
    if (matching_lsh(i) ~= matching_naive(i))
        wrong = wrong + 1;
    end
end
correct = nPoints - wrong;
correctRatio = correct/nPoints;

correctRatio





