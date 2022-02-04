% Simply create random orthogonal hyperplanes that sum to 1

% set seed
rng(0)

vectorLength = 128;

maxPlanes = 32;
hplanes = 2*rand(maxPlanes, vectorLength)-1;

for i = 1:maxPlanes
    
    % Other orthogonal vector that we can add to the hyperplane
    other = 2*rand(1, vectorLength) - 1;
   
    % Make it orthogonal to the hyperplanes
    for j = 1:(i-1)
        hplanes(i,:) = hplanes(i,:) ...
               - (hplanes(i,:) * hplanes(j,:)') * hplanes(j,:);
        
        other = other - (other * hplanes(j,:)') * hplanes(j,:);
    end
    
    % Make it sum to zero
    hplanes(i,:) = hplanes(i,:) ...
       - ((sum(hplanes(i,:)) / sum(other)) * other);
    
    % Normalize
    hplanes(i,:) = hplanes(i,:) / ...
        sqrt(hplanes(i,:) * hplanes(i,:)');
   
end


% Are they orthogonal?
cnt = 0;
for i = 1:32
    for j = (i+1):32
%         fprintf("hyperplanes(%d,:)*hyperplanes(%d,:)' = %f\n",...
%             i, j, hyperplanes(i,:)*hyperplanes(j,:)');
        if (abs(hplanes(i,:)*hplanes(j,:)') > 1e-10)
            cnt = cnt + 1;
        end
    end
end
fprintf("Number of hyperplanes not orthogonal: %d\n", cnt);

% What do the elements sum to?
cnt = 0;
for i = 1:32
%     fprintf("sum(hplanes(%d,:)) = %f\n", i, sum(hplanes(i,:)));
    if (abs(sum(hplanes(i,:))) > 1e-10)
        cnt = cnt + 1;
    end
end
fprintf("number of sums not zero: %d\n", cnt);


fileID = fopen('hyperplanes.dat','w');
fwrite(fileID,hplanes','float32');
fclose(fileID);

