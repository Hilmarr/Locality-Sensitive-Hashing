% Creating random hyperplanes that sum to 0

nPlanes = 1000;
vectorLength = 128;

hPlanes = 2*rand(2*nPlanes, vectorLength) - 1;

% for i = 1:nPlanes
%     sum1 = sum(hPlanes(i,:));
%     sqrdSize = hPlanes(i,:)*hPlanes(i,:)';
%     fprintf("%f, %f\n", sum1, sqrdSize);
% end

% Make them sum to zero
for i = 1:nPlanes
    sum1 = sum(hPlanes(i,:));
    sum2 = sum(hPlanes(i+nPlanes,:));
    hPlanes(i,:) = hPlanes(i,:) - (sum1/sum2) * hPlanes(i+nPlanes,:);
end

% Normalize vector
for i = 1:nPlanes
    hPlanes(i,:) = hPlanes(i,:) / sqrt(hPlanes(i,:)*hPlanes(i,:)');
end

hPlanes = hPlanes(1:1000,:);

% % Create normalized noise
% noise = 2*rand(nPlanes, vectorLength) - 1;
% for i = 1:nPlanes
%     noise(i,:) = noise(i,:) / sqrt(noise(i,:)*noise(i,:)');
% end
% 
% % Add noise and normalize
% for i = 1:nPlanes
%     hPlanes(i,:) = hPlanes(i,:) + 0.3 * noise(i,:);
%     hPlanes(i,:) = hPlanes(i,:) / sqrt(hPlanes(i,:)*hPlanes(i,:)');
% end

for i = 1:nPlanes
    sum1 = sum(hPlanes(i,:));
    sqrdSize = hPlanes(i,:)*hPlanes(i,:)';
    fprintf("%f, %f\n", sum1, sqrdSize);
end



fileID = fopen('hyperplanes.dat','w');
fwrite(fileID,hPlanes','float32');
fclose(fileID);

