function[area] = formsArea(forms,N)
%Determinung the Areas of each form
%forms - frame with the forms detected in the current frame
%N - number of forms detected in the current frame
area = zeros(N,1);
for j=1:N
    pos = find(forms == j);
    aux = zeros(size(forms));
    aux(pos) = 1;
    area(j) = bwarea(aux);
end
end