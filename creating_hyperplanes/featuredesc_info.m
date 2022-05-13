% Getting some basic info from the feature descriptors

% meanDistPoints = 536.3315;
% meanLengths = 508.6323;
% meanDistHyperplanes = 28.3056;
% meanDistPoints = 536.3315;
% Average distance between best matches = 187.751004


nPoints = 100000;
vectorLength = 128;

% Let these be the points to fit the hyperplanes to for now
% points = 2*rand(nPoints, vectorLength) - 1;
points = fvecs_read("../test_data/sift/sift_learn.fvecs", nPoints);

%% Find the average length of the feature descriptors
lengths = sqrt(sum(points.^2));
meanLengths = mean(lengths);

points = points';

%% Find average distance from hyperplanes (where elements sum to zero)
nPlanes = 50;
hyperplanes = rand(2*nPlanes, 128);
for i = 1:nPlanes
    hyperplanes(i,:) = hyperplanes(i,:) ...
            - ((sum(hyperplanes(i,:)) / sum(hyperplanes(nPlanes+i,:))) ...
              * hyperplanes(nPlanes+i,:));
end
hyperplanes = hyperplanes(1:nPlanes,:);
for i = 1:nPlanes
    hyperplanes(i,:) = hyperplanes(i,:) / sqrt(sum(hyperplanes(i,:).^2));
end

distAcc = 0;
for i = 1:nPlanes
    for j = 1:nPoints
        distAcc = distAcc + sqrt((points(j,:)*hyperplanes(i,:)')^2);
    end
end
meanDistHyperplanes = distAcc / (nPlanes*nPoints);

%% Find average distance between random points
distAcc = 0;
N = 200000;

for k = 1:N
    i = randi(nPoints,1);
    j = randi(nPoints,1);
    while j == i
        j = randi(nPoints,1);
    end
    distAcc = distAcc + sqrt(sum((points(i,:) - points(j,:)).^2));
end
meanDistPoints = distAcc / N;

    

% accDists = 0;
% 
% for i = 1:nPlanes
%     for j = 1:nPoints
%         accDists = accDists + points


