function keypoints = EdgeDetection(forms,N)
keypoints = cell(N,1);
fsize = size(forms);
for i=1:N
    pos = find(forms == i);
    aux = zeros(fsize);
    aux(pos) = 1;
    [~, threshold] = edge(aux, 'sobel');
    fudgeFactor = .5;
    edges = edge(aux,'sobel', threshold * fudgeFactor);
    [row,col] = find(edges == 1);
    keypoints{i} = [row,col];
end
end