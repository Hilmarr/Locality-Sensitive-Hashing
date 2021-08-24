
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

%% LSH
