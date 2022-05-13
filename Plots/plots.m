% data = [0.283, 0.009, 0.658, 1.01];
% names = ["hash base vectors", "create LSH tables"...
%     "Find potential matches", "match potential matches"];
% bar(data);
% set(gca,'xticklabel',names)
% set(gcf,'position',[100,100,768,512])
% ylabel("seconds")
% title("Run-time for each phase, loops optimized.")


% data = [129, 36.2, 33.9, 7.39, 6.65, 3.46];
% names = ["g++, CPU", "pgc++, CPU",...
%     "Parallelize hashing", "Parallelize matching", ...
%     "Loop optimization", "Memory movement in bulk"];
% bar(data);
% set(gca,'xticklabel',names)
% set(gcf,'position',[100,100,768,512])
% ylabel("seconds")
% title("Total run-time for each optimization.")
% 

% data = [61.9, 7.93, 4.19, 3.68, 0.726];
% names = ["g++, CPU",...
%     "pgc++, CPU",...
%     "Parallelize hashing", ...
%     "Loop optimization", ...
%     "Memory movement in bulk"];
% bar(data);
% set(gca,'xticklabel',names)
% set(gcf,'position',[100,100,768,512])
% ylabel("seconds")
% title("Hashing run-time for each optimization.")

data = [81.9, 33.5, 2.67, 2.44];
names = ["g++, CPU",...
    "pgc++, CPU",...
    "Parallelize matching", ...
    "Loop optimization"];
bar(data);
set(gca,'xticklabel',names)
set(gcf,'position',[100,100,700,512])
ylabel("seconds")
title("Total run-time for each optimization.")

% Serial using g++ 81.9s 2.25s 78.6s
% Serial using pgc++ 33.5s 0.257s 30.9s
% Parallelize 2.67s 0.304s 1.31s
% Loop optimization 2.44s 0.283s 1.01s

