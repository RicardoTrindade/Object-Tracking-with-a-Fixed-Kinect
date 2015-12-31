function cost_matrix = CostMatrixShape(forms,N,prev_forms,PN)
%forms - Matrix with the forms of the current frame
%N - number of forms in the current frame
%prev_forms - Matrix with the forms of the prevous frame
%PN - number of forms in the previous frame

cost_matrix = zeros(PN,N);
forms_aux = cell(N,1);
preforms_aux = cell(PN,1);
fsize = size(forms);

%Getting the frames and keyponts for each form and centering them
wind = 20;
for i=1:N
    pos = forms == i;
    aux = zeros(fsize);
    aux(pos) = 1;
    s = regionprops(aux,'centroid');
    c = s.Centroid;
    forms_aux{i} = aux(round(c(1,2)-wind):round(c(1,2)+wind),round(c(1,1)-wind):round(c(1,1)+wind));
end
for i=1:PN
    pos = prev_forms == i;
    aux = zeros(fsize);
    aux(pos) = 1;
    s = regionprops(aux,'centroid');
    c = s.Centroid;
    preforms_aux{i} = aux(round(c(1,2)-wind):round(c(1,2)+wind),round(c(1,1)-wind):round(c(1,1)+wind));
end

%Building the cost matrix
for i=1:N
    aux1 = forms_aux{i};
    kp_form = Edges(aux1);
    for j=1:PN
        aux2 = preforms_aux{j};
        kp_pform =  Edges(aux2);
        [index,match_distances] = matchFeatures(kp_form,kp_pform,'Unique',logical(true));
        cost = sum(match_distances)/size(index,1);%Normalizar de acordo com o numero de matches verificadas
        cost_matrix(j,i) = cost;
    end
end
cost_matrix = cost_matrix/(6*10^-6);
cost_matrix(isnan(cost_matrix))=0;
end