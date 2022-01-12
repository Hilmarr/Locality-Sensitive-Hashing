% Create orthogonal hyperplanes

vectorLength = 128;
nPlanes = 13;
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

for i = 1:nPlanes
    for j = (i+1):nPlanes
        fprintf("hyperplanes(%d,:)*hyperplanes(%d,:)' = %f\n",...
            i, j, hyperplanes(i,:)*hyperplanes(j,:)');
    end
end