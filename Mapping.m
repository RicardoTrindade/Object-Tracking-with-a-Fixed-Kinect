function map = Mapping(tracks,frame,k,c)

stop = 1;
s = length(tracks);
if k == 0
aux2 = tracks(frame).undetect(c);
map(frame) = aux2;
j = frame+1;
else
    aux2 = tracks(frame).match(1,2);
    map(frame) = aux2;
    j = frame+1;
end
for i=j:s
    aux = tracks(i).match;
    if isempty(aux) && stop==0
        continue;
    end
    if (isempty(aux) || isempty(find(aux(:,1)==aux2, 1))) && stop==1
%         if isempty(find(tracks(i).undetect(:)==aux2, 1)) 
        break;
%         else
%             map(i) = aux2;
%         end
    else
    if ~isempty(aux)
        stop = 1;
        [lin,~] = find(aux(:,1) == aux2);
        aux2 = aux(lin,2);
%         disp(i);
        map(i) = aux2;
    end
    end
end
    
end
