function cost_matrix = CostMatrixAreas(area,prev_area,N)
%Determining the Cost Matrix for the areas
%area - area of the forms in the current frame
%prev_area - area of the forms in the previous frame
%N - number of forms in the current frame
cost_matrix = zeros(size(prev_area,1),N);
for i = 1:N
    for j1 = 1:size(prev_area,1)
        error = (area(i)-prev_area(j1))^2;
        cost_matrix(j1,i) = error;
    end
end
cost_matrix = cost_matrix/196;
cost_matrix(isnan(cost_matrix)) = 0;
end