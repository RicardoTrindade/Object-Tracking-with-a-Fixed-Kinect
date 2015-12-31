function [nframe,kappas] = BeginTrack(tracks)

k = 0;
s = length(tracks);
nframe = [];
kappas = [];
j = 1;
aux2 = tracks(j).match;
if ~isempty(aux2)
    k = 1;
    for a = size(aux2,1)
        nframe(end+1) = j;
        kappas(end+1) = k;
    end
else
aux1 = tracks(j).undetect;
while j<=s && isempty(aux1) 
    j = j+1;
    aux1 = tracks(j).undetect; 
end
nframe(end+1) = j;
kappas(end+1) = k;
end
k = 0;
j = j+1;
for i=j:s
    aux = tracks(i).match;
%     if~isempty(aux) && ~all(aux(:,1)==aux(:,2),1) && ~isempty(tracks(i).undetect)
        if~isempty(aux) && ~isempty(tracks(i).undetect)
        leng = size(tracks(i).undetect,1);
        if leng>1
            for i1=1:leng
                nframe(end+1) = i;
                kappas(end+1) = k;
            end
        else
        nframe(end+1) = i;
        kappas(end+1) = k;
        end
    end
end
end