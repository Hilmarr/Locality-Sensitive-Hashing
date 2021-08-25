
%% set up the points
clear all; close all; clc;

nPoints = 100;
vectorLength = 128;
noiseScale = 0.2;
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

matching_naive = zeros(nPoints, 2);

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
    matching_naive(i,:) = [i, match];
end


%% Contstruct LSH table (one LSH table)

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

% % temporary storage as a number of linked lists
% % idx: point, next
% % next == 0, represents that it is a tail node,
% % point == 0, represents that there is no point stored there
% % This works in matlab since indexing starts from 1, in other languages,
% % these metavalues should be -1 or some other invalid index.
% groupTable = zeros(nPoints, 2);

% array of points grouped by what box they got hashed into
groupArray = zeros(nPoints, 1);

% Put points1 into our lsh table
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
    
    
% LATER: check for duplicates
% (in order to check whether using hyperplanes that go through the origin
%  is better than more random hyperplanes, should probably plot number
%  of duplicates given the length of the hyperplanes)

%% Match points2 with points1 using LSH hash table
















