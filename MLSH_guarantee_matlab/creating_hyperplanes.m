% set seed
rng(0)

%% Create orthogonal hyperplanes

vectorLength = 128;
nPlanes = 128;
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

% % Are they orthogonal? Yes
% for i = 1:nPlanes
%     for j = (i+1):nPlanes
%         fprintf("hyperplanes(%d,:)*hyperplanes(%d,:)' = %f\n",...
%             i, j, hyperplanes(i,:)*hyperplanes(j,:)');
%     end
% end

%% Normalize

for i = 1:13
    h = hyperplanes(i,:);
    hyperplanes(i,:) = h / sqrt(h*h');
end

% nNewPlanes = nPlanes - 13;
% newPlanes = rand(nNewPlanes, vectorLength);

%%  Calculate random hyperplanes that are also normal to the previous ones

for i = 14:nPlanes
    new = rand(1, vectorLength);
    for j = 1:(i-1)
        new = new - (new * hyperplanes(j,:)') * hyperplanes(j,:);
    end
    hyperplanes(i,:) = new / sqrt(new*new');
end

% % Are they orthogonal? Yes
% for i = 1:nPlanes
%     for j = (i+1):nPlanes
%         fprintf("hyperplanes(%d,:)*hyperplanes(%d,:)' = %f\n",...
%             i, j, hyperplanes(i,:)*hyperplanes(j,:)');
%     end
% end
% 
% % What do the elements sum to?
% for i = 1:nPlanes
%     fprintf("sum(hyperplanes(%d,:)) = %f\n", i, sum(hyperplanes(i,:)));
% end

%% Find the ones with equal sum distribution

nNewPlanes = nPlanes - 13;
inner_sum = zeros(nNewPlanes,1);
abs_inner_sum = zeros(nNewPlanes,1);

for i = 1:nNewPlanes
    inner_sum(i) = sum(hyperplanes(i+13,:));
    abs_inner_sum(i) = abs(inner_sum(i));
end

[abs_inner_sum2, index] = sort(abs_inner_sum);

inner_sum2 = zeros(nNewPlanes,1);
new_hyperplanes = zeros(nNewPlanes,vectorLength);

for i = 1:nNewPlanes
    new_hyperplanes(i,:) = hyperplanes(index(i)+13,:);
    inner_sum2(i) = inner_sum(index(i));
end

% % Are they orthogonal? Yes
% cnt = 0;
% for i = 1:57
%     for j = (i+1):57
% %         fprintf("new_hyperplanes(%d,:)*new_hyperplanes(%d,:)' = %f\n",...
% %             i, j, new_hyperplanes(i,:)*new_hyperplanes(j,:)');
%         if (abs(new_hyperplanes(i,:)*new_hyperplanes(j,:)') > 1e-10)
%             cnt = cnt + 1;
%         end
%     end
% end
% fprintf("new_hyperplanes not orthogonal: %d\n", cnt);

% % Correct sorting? Yes
% fprintf("%f\n", max(abs(sum(new_hyperplanes(i,:)) - inner_sum2(i))));

%% Make the hyperplanes sum to 0 while still being orthogonal
%  to the original hyperplanes

for i = 1:57
    new_hyperplanes(i,:) = new_hyperplanes(i,:) ...
        - ((inner_sum2(i) / inner_sum2(57+i))...
          * new_hyperplanes(57+i,:));
end

% % Are they orthogonal? Yes
% cnt = 0;
% for i = 1:57
%     for j = (i+1):57
% %         fprintf("new_hyperplanes(%d,:)*new_hyperplanes(%d,:)' = %f\n",...
% %             i, j, new_hyperplanes(i,:)*new_hyperplanes(j,:)');
%         if (abs(new_hyperplanes(i,:)*new_hyperplanes(j,:)') > 1e-10)
%             cnt = cnt + 1;
%         end
%     end
% end
% fprintf("new_hyperplanes not orthogonal: %d\n", cnt);

inner_sum3 = zeros(nNewPlanes,1);
for i = 1:nNewPlanes
    inner_sum3(i) = sum(new_hyperplanes(i,:));
end

new_hyperplanes2 = new_hyperplanes(1:57,:);

% Normalize again
for i = 1:57
    h = new_hyperplanes2(i,:);
    new_hyperplanes2(i,:) = h / sqrt(h*h');
end

positiveRatio = zeros(57,1);
for i = 1:57
    positiveRatio(i) = sum(new_hyperplanes2(i,:) > 0) / 128;
end
positiveRatio = positiveRatio - mean(positiveRatio);
positiveRatio = positiveRatio / std(positiveRatio);
positiveRatio = abs(positiveRatio);

my_medians = zeros(57,1);
for i = 1:57
    my_medians(i) = ...
        abs(median(new_hyperplanes2(i,:)));
end
my_medians = my_medians - mean(my_medians);
my_medians = my_medians / std(my_medians);
my_medians = abs(my_medians);

% Create some cost function that I want to sort by, I only really need some
% of them. Just want the distribution to be compact if possible, not really
% too important.
cost = positiveRatio + my_medians;

[r, index] = sort(cost);

new_hyperplanes3 = zeros(57, 128);

for i = 1:57
    new_hyperplanes3(i,:) = new_hyperplanes2(index(i),:);
end

% % Are they orthogonal? Yes
% cnt = 0;
% for i = 1:57
%     for j = (i+1):57
% %         fprintf("new_hyperplanes2(%d,:)*new_hyperplanes2(%d,:)' = %f\n",...
% %             i, j, new_hyperplanes2(i,:)*new_hyperplanes2(j,:)');
%         if (abs(new_hyperplanes3(i,:)*new_hyperplanes3(j,:)') > 1e-10)
%             cnt = cnt + 1;
%         end
%     end
% end
% fprintf("new_hyperplanes3 not orthogonal: %d\n", cnt);



%% Add them to the final hyperplanes

final_hplanes = zeros(70,128);
final_hplanes(1:13,:) = hyperplanes(1:13,:);
final_hplanes(14:70,:) = new_hyperplanes3(:,:);

% Are they orthogonal? Yes
cnt = 0;
for i = 1:70
    for j = (i+1):70
%         fprintf("final_hplanes(%d,:)*final_hplanes(%d,:)' = %f\n",...
%             i, j, final_hplanes(i,:)*final_hplanes(j,:)');
        if (abs(final_hplanes(i,:)*final_hplanes(j,:)') > 1e-10)
            cnt = cnt + 1;
        end
    end
end
fprintf("Number of hyperplanes not orthogonal: %d\n", cnt);

% What do the elements sum to?
cnt = 0;
for i = 1:70
%     fprintf("sum(final_hplanes(%d,:)) = %f\n", i, sum(final_hplanes(i,:)));
    if (abs(sum(final_hplanes(i,:))) > 1e-10)
        cnt = cnt + 1;
    end
end
fprintf("number of sums not zero: %d\n", cnt);


fileID = fopen('hyperplanes.dat','w');
fwrite(fileID,final_hyperplanes,'float32');
fclose(fileID);