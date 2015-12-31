function keypoints = Edges(form)
%Determine the points of the Edges for the form

[~, threshold] = edge(form, 'sobel');
fudgeFactor = .5;
edges = edge(form,'sobel', threshold * fudgeFactor);
[row,col] = find(edges == 1);
keypoints = [row,col];
end