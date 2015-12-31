function [filtforms,newN] = AreaFilt(forms,N)
%Filter for remotion of very small structures
%forms - frame with all the forms
%N - number of forms in the frame
filtforms = zeros(size(forms));
threshold = 7;%Minimum area for the form to be accepted
areas = formsArea(forms,N);
newN = 0;
j = 1;

for i=1:N
    if areas(i)<=threshold
    else
        pos = forms == i;
        area = zeros(size(forms));
        area(pos) = j;
        filtforms = filtforms+area;
        j = j+1;
        newN = newN+1;
    end
end
end