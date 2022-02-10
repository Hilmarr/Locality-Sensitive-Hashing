% set seed
rng(0)

nPoints = 25000;
vectorLength = 128;

% Let these be the points to fit the hyperplanes to for now
% points = 2*rand(nPoints, vectorLength) - 1;
points = fvecs_read("../test_data/sift/sift_learn.fvecs", [nPoints, 2*nPoints]);
points = points';

maxPlanes = 20;
hyperplanes = zeros(maxPlanes, vectorLength);

nAlternatives = 100;

splitScore = zeros(1, nAlternatives);
absDists = zeros(1,nAlternatives);

groupSizes = zeros(1,1);
groupSizes(1) = nPoints;
groupArray = 1:nPoints;
groupArray2 = zeros(1,nPoints); % Temporary array
whichSide = zeros(1,nPoints); % Temporary array

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
    nGroups = 2^(nPlanes-1);
    
    for i = 1:nAlternatives
        splitScore(i) = 0;
        idx = 1;
        for j = 1:nGroups
            positiveSide = 0;
            for k = idx:(idx+groupSizes(j)-1)
                prod = alts(i,:) * points(groupArray(k),:)';
                if (prod > 0)
                    positiveSide = positiveSide + 1;
                end
                absDists(i) = absDists(i) + abs(prod);
            end
            negativeSide = groupSizes(j) - positiveSide;
            splitScore(i) = splitScore(i) + min(negativeSide, positiveSide);
            
            idx = idx + groupSizes(j);
        end
    end

    splitScore = splitScore / max(splitScore);
    
    absDists = absDists / max(absDists);

    score = absDists + splitScore;
    
    % Add hyperplane with best scores to the hyperplane list
    [minimum_score, bestIdx] = max(score);
    
    hyperplanes(nPlanes,:) = alts(bestIdx,:);
    
    % Recreate groupSizes and groupArray
    hplane = hyperplanes(nPlanes,:);
    nGroups2 = 2*nGroups;
    groupSizes2 = zeros(1,nGroups2);
    
    idx = 1;
    for i = 1:nGroups
        for k = idx:(idx+groupSizes(i)-1)
            prod = hplane * points(groupArray(k),:)';
%             fprintf("prod=%f, groupArray(%d)=%d,  ", prod, i, groupArray(i));
            if (prod <= 0)
%                 fprintf("negative side\n");
                groupSizes2(2*i-1) = groupSizes2(2*i-1) + 1;
                whichSide(k) = 1;
            else
%                 fprintf("positive side\n");
                groupSizes2(2*i) = groupSizes2(2*i) + 1;
                whichSide(k) = 2;
            end
        end
        nextIdx1 = idx;
        nextIdx2 = idx + groupSizes2(2*i-1);
        for k = idx:(idx+groupSizes(i)-1)
            if (whichSide(k) == 1)
                groupArray2(nextIdx1) = groupArray(k);
                nextIdx1 = nextIdx1 + 1;
            else
                groupArray2(nextIdx2) = groupArray(k);
                nextIdx2 = nextIdx2 + 1;
            end
        end
        
        idx = idx + groupSizes(i);
%         fprintf("idx=%d\n", idx);
    end
    
%     for i = 1:nPoints
%         groupArray(i) = groupArray2(i);
%     end

    groupArray = groupArray2;
    
    groupSizes = groupSizes2;
end

fprintf("Finished finding hyperplanes\n");

% Test

positiveSide = zeros(1,maxPlanes);
absDists = zeros(1,maxPlanes);

positiveSide2 = zeros(1,maxPlanes);
squaredDists2 = zeros(1,maxPlanes);

% Random hyperplanes
% Need to make them sum to one when testing on real points
randplanes = rand(maxPlanes*2, vectorLength);

% Normalize random hyperplanes
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
        absDists(i) = absDists(i) + prod*prod;
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
absDists = absDists / nPoints;
dists = sqrt(absDists);
positiveSide2 = positiveSide2 / nPoints;
squaredDists2 = squaredDists2 / nPoints;
dists2 = sqrt(squaredDists2);

% % Are they orthogonal?
% cnt = 0;
% for i = 1:maxPlanes
%     for j = (i+1):maxPlanes
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

