% set seed
rng(0)

nPoints = 25000;
vectorLength = 128;

% Let these be the points to fit the hyperplanes to for now
% points = 2*rand(nPoints, vectorLength) - 1;
points = fvecs_read("../test_data/sift/sift_learn.fvecs", nPoints);
points = points';

maxPlanes = 32;
hyperplanes = zeros(maxPlanes, vectorLength);

nAlternatives = 1000;

positiveSide = zeros(1,nAlternatives);
squaredDists = zeros(1,nAlternatives);

for nPlanes = 1:maxPlanes
    fprintf("Finding hyperplane %d of %d\n", nPlanes, maxPlanes);
    
    % Get alternatives as well as some extra to adjust them
    alts = 2*rand(2*nAlternatives, vectorLength) - 1;
   
    % Make them orthogonal to the hyperplanes
    for i = 1:(2*nAlternatives)
        for j = 1:(nPlanes-1)
            alts(i,:) = alts(i,:) - ...
                (alts(i,:) * hyperplanes(j,:)') * hyperplanes(j,:);
        end
    end
    
    % Make them all sum to zero
    % (Only a good idea if the points are all positive as is the case
    %  for the actual sift descriptors later, would not be necessary if one
    %  could try even more hyperplanes, but might shorten the search)
    for i = 1:nAlternatives
        alts(i,:) = alts(i,:) ...
            - ((sum(alts(i,:)) / sum(alts(nAlternatives+i,:))) ...
              * alts(nAlternatives+i,:));
    end
    
    % Normalize vectors
    for i = 1:nAlternatives
        alts(i,:) = alts(i,:) / sqrt(alts(i,:) * alts(i,:)');
    end
   
    % Create scores for how well it separates the search space
    % and for average distance between point and hyperplane
    for i = 1:nAlternatives
        positiveSide(i) = 0;
        squaredDists(i) = 0;
        for j = 1:nPoints
            prod = alts(i,:) * points(j,:)';
            if (prod > 0)
                positiveSide(i) = positiveSide(i) + 1;
            end
            squaredDists(i) = squaredDists(i) + prod*prod;
        end
    end
    positiveSide = positiveSide / nPoints;
    positiveSide = (0.5 - positiveSide).^2;
    positiveSide = positiveSide / std(positiveSide);
    
    squaredDists = squaredDists / nPoints;
    for i = 1:nAlternatives
        squaredDists(i) = 1 / (1 + squaredDists(i));
    end
    squaredDists = squaredDists / std(squaredDists);

    cost = positiveSide + squaredDists;
    
    % Add hyperplane with best scores to the hyperplane list
    [minimum_score, bestIdx] = min(cost);
    
    hyperplanes(nPlanes,:) = alts(bestIdx,:);
    
end

fprintf("Finished finding hyperplanes\n");

% Test

positiveSide = zeros(1,maxPlanes);
squaredDists = zeros(1,maxPlanes);

positiveSide2 = zeros(1,maxPlanes);
squaredDists2 = zeros(1,maxPlanes);

% Random hyperplanes
% Need to make them sum to one when testing on real points
randplanes = rand(maxPlanes*2, vectorLength);

% Make random hyperplanes sum to zero
for i = 1:maxPlanes
    randplanes(i,:) = randplanes(i,:) ...
        - (sum(randplanes(i,:)) / sum(randplanes(maxPlanes+i,:))) ...
          * randplanes(maxPlanes+i,:);
end

% Normalize random hyperplanes
for i = 1:maxPlanes
    randplanes(i,:) = randplanes(i,:) / ...
        sqrt(randplanes(i,:) * randplanes(i,:)');
end


% Create scores for how well it separates the search space
% and for average distance between point and hyperplane
for i = 1:maxPlanes
    for j = 1:nPoints
        prod = hyperplanes(i,:) * points(j,:)';
        if (prod > 0)
            positiveSide(i) = positiveSide(i) + 1;
        end
        squaredDists(i) = squaredDists(i) + prod*prod;
    end
end
    
for i = 1:maxPlanes
    for j = 1:nPoints
        prod = randplanes(i,:) * points(j,:)';
        if (prod > 0)
            positiveSide2(i) = positiveSide2(i) + 1;
        end
        squaredDists2(i) = squaredDists2(i) + prod*prod;
    end
end

positiveSide = positiveSide / nPoints;
squaredDists = squaredDists / nPoints;
dists = sqrt(squaredDists);
positiveSide2 = positiveSide2 / nPoints;
squaredDists2 = squaredDists2 / nPoints;
dists2 = sqrt(squaredDists2);

% % Are they orthogonal?
% cnt = 0;
% for i = 1:32
%     for j = (i+1):32
% %         fprintf("hyperplanes(%d,:)*hyperplanes(%d,:)' = %f\n",...
% %             i, j, hyperplanes(i,:)*hyperplanes(j,:)');
%         if (abs(hyperplanes(i,:)*hyperplanes(j,:)') > 1e-10)
%             cnt = cnt + 1;
%         end
%     end
% end
% fprintf("Number of hyperplanes not orthogonal: %d\n", cnt);

fprintf("DONE\n");

fileID = fopen('hyperplanes.dat','w');
fwrite(fileID,hyperplanes','float32');
fclose(fileID);
% 
% 
% fid = fopen('hyperplanes.dat','r');
% hyperplanes = fread(fid, [128, 32], 'float32');
% hyperplanes = hyperplanes';
% fclose(fid)

