
%% set up the points
clear all; close all; clc;

nPoints = 5000;
vectorLength = 128;
noiseScale = 0.5;
points1 = 2*rand(nPoints, vectorLength) - 1;
% Noise generated from a uniform random distribution.
noise = (2*rand(nPoints, vectorLength) - 1) * noiseScale;
% Noise generated from a normal distribution
% noise = randn(nPoints, vectorLength) * noiseScale;
points2 = points1 + noise;

% Normalize the vectors so their absolute values are all 1 (l2 norms are 1)
for i = 1:nPoints
    points1(i,:) = points1(i,:) / norm(points1(i,:), 2);
    points2(i,:) = points2(i,:) / norm(points2(i,:), 2);
end

%% LSH
lshStart = tic;

nTables = 3;
nPlanes = round(log2(nPoints));

% Store best matches for LSH
matching_lsh = zeros(nPoints, 1);
% Store the distance for each best match
bestMatchDists = zeros(nPoints) + 1e10;
% Booleans
%   True if the point was within some tolerance of hyperplane
%   False otherwise
closeToHP = zeros(nPoints, nPlanes);
tol = 1e-1;
tol = tol*tol;

% Initialize LSH tables
nBoxes = 2^nPlanes;
indexGroupMap = zeros(nPoints, 1);
groupSizeMap = zeros(nBoxes, 1);
groupIndexMap = zeros(nBoxes, 1);
groupIndexMapTails = zeros(nBoxes, 1);

for table = 1:nTables
    
    % -- Construct LSH table --

    % Initialize random hyperplanes
    hyperplanes = 2*rand(nPlanes, vectorLength) - 1;
    
    % Normalize hyperplanes
    for i = 1:nPlanes
        hyperplanes(i,:) = hyperplanes(i,:) / norm(hyperplanes(i,:), 2);
    end

    % Initialize LSH tables
    indexGroupMap = 0 * indexGroupMap;
    groupSizeMap = 0 * groupSizeMap;
    groupIndexMap = 0 * groupIndexMap;
    groupIndexMapTails = 0 * groupIndexMapTails;

    % array of points grouped by what box they got hashed into
    groupArray = zeros(nPoints, 1);

    % - Put points1 into our lsh table -

    % Calculate hash values, keep track of sizes of each group
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


    % -- Match points2 with points1 using LSH hash table --

    % Calculating hash codes happens separately from searching,
    % this makes the program easier to parallelize at some later point
    
    closeToHP = 0 * closeToHP;

    % Calculate hash values
    hashValuesStart2 = tic;
    for i = 1:nPoints
        point = points2(i,:);
        hashcode = 0;
        for j = 1:nPlanes
            hplane = hyperplanes(j,:)';
            tmp = point*hplane;
            if (tmp > 0)
                hashcode = bitor(hashcode, bitshift(1, j-1));
            end
            closeToHP(i, j) = tmp < tol;
%             if (tmp*tmp < tol) 
%                 fprintf("tmp*tmp=%f, tol=%f\n", tmp*tmp, tol);
%             end
        end
        hashcode = hashcode+1; % Because matlab is weird with indexing
        indexGroupMap(i) = hashcode;
    end

    for i = 1:nPoints
        % Find table and elements in points1 to match with
        changed = 0;
        bestMatchDist = bestMatchDists(i);
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
            if (diff < bestMatchDist)
                bestMatchDist = diff;
                match = idx;
                changed = 1;
            end
        end
        
        for k = 1:nPlanes
            if (closeToHP(i, k))
                hashcode2 = bitxor(hashcode-1, bitshift(1, k-1))+1;

                size = groupSizeMap(hashcode2);
                startIdx = groupIndexMap(hashcode2);

                for j = startIdx:(startIdx+size-1)
                    idx = groupArray(j);
                    diff = sum((points2(i,:) - points1(idx,:)) .^ 2);
                    if (diff < bestMatchDist)
                        bestMatchDist = diff;
                        match = idx;
                        changed = 1;
                    end
                end
            end
        end
        
        if (changed == 1)
            matching_lsh(i) =  match;
            bestMatchDists(i) = bestMatchDist;
        end
    end

end

lshTime = toc(lshStart);


%% Deal with time

wrong = 0;
for i = 1:nPoints
    if (matching_lsh(i) ~= i)
        wrong = wrong + 1;
    end
end
correct = nPoints - wrong;
correctRatio = correct/nPoints;


fprintf("N = %d\n", nPoints);
fprintf("Total matching time using LSH:         %f\n", lshTime);
fprintf("LSH correct ratio:                  %f\n", correctRatio);