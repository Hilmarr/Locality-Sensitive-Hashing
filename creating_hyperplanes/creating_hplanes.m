% set seed
rng(0)

nPoints = 100;
vectorLength = 128;

% Let these be the points to fit the hyperplanes to for now
points = 2*rand(nPoints, vectorLength) - 1;

maxPlanes = 32;
hyperplanes = zeros(maxPlanes, vectorLength);

nAlternatives = 100;

positiveSide = zeros(1,nAlternatives);
squaredDists = zeros(1,nAlternatives);

for nPlanes = 1:maxPlanes
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
    
    for i = 1:nAlternatives
        squaredDists(i) = nPoints / (1 + squaredDists(i)*squaredDists(i));
    end
    squaredDists = squaredDists / std(squaredDists);

    cost = positiveSide + squaredDists;
    
    % Add hyperplane with best scores to the hyperplane list
    [minimum_score, bestIdx] = min(cost);
    
    hyperplanes(nPlanes,:) = alts(bestIdx,:);
    
end


% Test

positiveSide = zeros(1,maxPlanes);
squaredDists = zeros(1,maxPlanes);

positiveSide2 = zeros(1,maxPlanes);
squaredDists2 = zeros(1,maxPlanes);

% Random hyperplanes
randplanes = rand(maxPlanes, vectorLength);

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
positiveSide2 = positiveSide2 / nPoints;
squaredDists2 = squaredDists2 / nPoints;

% Are they orthogonal?
cnt = 0;
for i = 1:32
    for j = (i+1):32
%         fprintf("hyperplanes(%d,:)*hyperplanes(%d,:)' = %f\n",...
%             i, j, hyperplanes(i,:)*hyperplanes(j,:)');
        if (abs(hyperplanes(i,:)*hyperplanes(j,:)') > 1e-10)
            cnt = cnt + 1;
        end
    end
end
fprintf("Number of hyperplanes not orthogonal: %d\n", cnt);




