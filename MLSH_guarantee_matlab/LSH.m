clear all; close all; clc;

% set seed
rng(0)

%% Create orthogonal hyperplanes

vectorLength = 128;
nPlanes = 13;
hyperplanes = zeros(nPlanes, vectorLength);

for i = 1:7
    val = 1;
    steps = 2^i;
    stepLen = vectorLength / steps;
    for j = 1:stepLen:vectorLength
        for k = j:(j+stepLen-1)
            hyperplanes(i,k) = val;
        end
        val = -val;
    end
end

for i = 1:6
    for j = 0:(vectorLength-1)
        hyperplanes(i+7, j+1) = ...
        hyperplanes(i, mod(j+(2^(6-i)), vectorLength) + 1);
    end
end

% for i = 1:nPlanes
%     for j = (i+1):nPlanes
%         fprintf("hyperplanes(%d,:)*hyperplanes(%d,:)' = %f\n",...
%             i, j, hyperplanes(i,:)*hyperplanes(j,:)');
%     end
% end

for i = 1:nPlanes
    h = hyperplanes(i,:);
    hyperplanes(i,:) = h / sqrt(h*h');
end


%% set up the points

nPlanes = 8;
nPoints = 100;
vectorLength = 128;
noiseScale = 0.1;
points1 = 2*rand(nPoints, vectorLength) - 1;
% Noise generated from a uniform random distribution.
noise = (2*rand(nPoints, vectorLength) - 1) * noiseScale;
% Noise generated from a normal distribution
% noise = randn(nPoints, vectorLength) * noiseScale;
points2 = points1 + noise;


%% Construct LSH table (one LSH table)
lshStart = tic;

% % Give the hyperplanes the length of somewhere between 0 and maxHPLength
% % Not needed if all hyperplanes go through the origin
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
        if (point*hplane > 0)
            hashcode = bitor(hashcode, bitshift(1, j-1));
        end
    end
    hashcode = hashcode+1; % Because matlab is weird with indexing
    indexGroupMap(i) = hashcode;
    groupSizeMap(hashcode) = groupSizeMap(hashcode)+1;
end
hashValuesTime1 = toc(hashValuesStart1);

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
    
    % Increment the index tail mapping
    groupIndexMapTails(hashcode) = idx + 1;
    
    groupArray(idx) = i;
    
end

tablesCreatedTime = toc(lshStart);


%% Match points2 with points1 using LSH hash table
matchingStart = tic;

threshold = 0.1;
threshold = threshold*threshold; % Squared max distance
% Squared dists from hyperplanes
sqrdDists = zeros(nPoints, nPlanes);

% Calculating hash codes happens separately from searching,
% this makes the program easier to parallelize at some later point

% Calculate hash values
hashValuesStart2 = tic;
for i = 1:nPoints
    point = points2(i,:);
    hashcode = 0;
    for j = 1:nPlanes
        hplane = hyperplanes(j,:)';
        proj = point*hplane;
        if (proj > 0)
            hashcode = bitor(hashcode, bitshift(1, j-1));
        end
        sqrdDists(i,j) = proj*proj;
    end
    hashcode = hashcode+1; % Because matlab is weird with indexing
    indexGroupMap(i) = hashcode;
end
hashValuesTime2 = toc(hashValuesStart2);

% for i = 1:nPoints
%     fprintf("Point %d:\n", i);
%     for j = 1:nPlanes
%         fprintf("  dist from h%d = %f\n", j, sqrt(sqrdDists(i,j)));
%     end
%     fprintf("\n");
% end

% Should later be naturally extended like an arraylist
groupMapExt = zeros(nPoints*5,1);
groupMapExtIndices = zeros(nPoints + 1);

cnt = 1;
for i = 1:nPoints
    hashcode = indexGroupMap(i) - 1;
%     fprintf("hashcode%d = %d\n", i, hashcode);
    groupMapExtIndices(i) = cnt;
    groupMapExt(cnt) = hashcode + 1;
    cnt = cnt + 1;
    % Replace this with a recursive function that can also find
    % combination of these values
    for j = 1:nPlanes
        if (sqrdDists(i,j) < threshold)
            groupMapExt(cnt) = bitxor(hashcode, bitshift(1, j-1)) + 1;
%             fprintf("  j=%d, hashcode=%d\n", j, groupMapExt(cnt));
            cnt = cnt + 1;
        end
    end
%     fprintf("\n");
end
groupMapExtIndices(nPoints + 1) = cnt;

% for i = 1:nPoints
%     fprintf("Point %d hash codes:\n", i);
% %     for j = groupMapExtIndices(i):(groupMapExtIndices(i+1)-1)
% %         fprintf("  %s\n", dec2bin(groupMapExt(j), nPlanes));
% %     end
%     for j = groupMapExtIndices(i):(groupMapExtIndices(i+1)-1)
%         fprintf("  %d\n", groupMapExt(j));
%     end
%     fprintf("\n");
% end



% Match the points using lsh
matching_lsh = zeros(nPoints, 1);

for i = 1:nPoints
    
    shortestDist = 1e10;
    
    for k = groupMapExtIndices(i):(groupMapExtIndices(i+1)-1)
        
        % Find table and elements in points1 to match with
        hashcode = groupMapExt(k);
        size = groupSizeMap(hashcode);
        startIdx = groupIndexMap(hashcode);

        % Match the points
        for j = startIdx:(startIdx+size-1)
            % Take euclidean distance between the vector
            idx = groupArray(j);
            diff = sum((points2(i,:) - points1(idx,:)) .^ 2);
            % If euclidean distance is better than the best match,
            % we have a new best match
            if (diff < shortestDist)
                shortestDist = diff;
                match = idx;
%                 fprintf("i=%d, match=%d\n",i,idx);
            end
        end
    end
    
    if (shortestDist ~= 1e10)
        matching_lsh(i) =  match;
    end
end

matchingTime = toc(matchingStart);
lshTime = toc(lshStart);
hashValuesTime = hashValuesTime1 + hashValuesTime2;



%% Deal with time

% naiveMatchingTime
% tablesCreatedTime
% matchingTime
% lshTime
% hashValuesTime1
% hashValuesTime2
% hashValuesTime


wrong = 0;
for i = 1:nPoints
    if (matching_lsh(i) ~= i)
        wrong = wrong + 1;
    end
end
correct = nPoints - wrong;
correctRatio = correct/nPoints;

% correctRatio

fprintf("N = %d\n", nPoints);
% fprintf("Matching time using naive method:      %f\n", naiveMatchingTime);
fprintf("Total matching time using LSH:         %f\n", lshTime);
disp(" *");
fprintf("  time portion: hyperplane,vector multiplication:  %f%% \n", ...
    100*hashValuesTime/lshTime);
disp(" *");
fprintf("  time portion spent creating and filling tables:  %f%% \n", ...
    100*tablesCreatedTime/lshTime);
fprintf("  time portion spent matching 2nd dataset:         %f%% \n", ...
    100*matchingTime/lshTime);
fprintf("\n\n");
% fprintf("Time ratio (lsh time / naive time): %f\n", (lshTime/naiveMatchingTime));
fprintf("LSH correct ratio:                  %f\n", correctRatio);


