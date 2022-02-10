nPlanes = 8;

fileID = fopen('hyperplanes.dat','r');
hplanes = fread(fileID,[128,nPlanes],'float32');
fclose(fileID);
hplanes = hplanes';

% Are they orthogonal?
cnt = 0;
for i = 1:nPlanes
    for j = (i+1):nPlanes
%         fprintf("hyperplanes(%d,:)*hyperplanes(%d,:)' = %f\n",...
%             i, j, hyperplanes(i,:)*hyperplanes(j,:)');
        if (abs(hplanes(i,:)*hplanes(j,:)') > 1e-6)
            cnt = cnt + 1;
        end
    end
end
fprintf("Number of hyperplanes not orthogonal: %d\n", cnt);

% What do the elements sum to?
cnt = 0;
for i = 1:nPlanes
%     fprintf("sum(hplanes(%d,:)) = %f\n", i, sum(hplanes(i,:)));
    if (abs(sum(hplanes(i,:))) > 1e-6)
        cnt = cnt + 1;
    end
end
fprintf("number of sums not zero: %d\n", cnt);

nPoints = 25000;
points = fvecs_read("../test_data/sift/sift_learn.fvecs", nPoints);
points = points';

nBoxes = 2^nPlanes;
groupSizeMap = zeros(1, nBoxes);
% hplaneDists = zeros(nPoints, nPlanes);
distSum = 0;

for i = 1:nPoints
    point = points(i,:);
    hashcode = 0;
    for j = 1:nPlanes
        hplane = hplanes(j,:)';
        dist = point*hplane;
        if (dist > 0)
            hashcode = bitor(hashcode, bitshift(1, j-1));
        end
        distSum = distSum + abs(dist);
%         hplaneDists(i, j) = dist;
    end
    hashcode = hashcode+1; % Because matlab is weird with indexing
    groupSizeMap(hashcode) = groupSizeMap(hashcode)+1;
end



avgDist = distSum / (nPoints * nPlanes);

fprintf("avgDist = %f\n", avgDist);
fprintf("std(groupSizeMap) = %f\n", std(groupSizeMap));
fprintf("max(groupSizeMap) = %f\n", max(groupSizeMap));