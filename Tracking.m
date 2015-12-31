clear
close all
clc
load('tracks.mat');


n=length(tracks);
caux_struct=tracks;
obj_count=1;
for i=38:50
    
    if isempty(aux_struct(i).match)
        continue;
    end
    %     if isempty(aux_struct(i+1).match)
    %         continue;
    %     end
    for j=1:size(aux_struct(i).match(:,2),1)
        prev=aux_struct(i).match(j,2);
        k=i+1;
        [gone]=find(aux_struct(k).unassigned==prev);
        if(~isempty(gone) || isempty(aux_struct(k).match))
            obj_count=obj_count+1;
            break;
        end
        next=aux_struct(k).match(:,1);
        tracker=find(next==prev);
        if(~isempty(tracker))
            aux_struct(i).match(j,1:2)
            aux_struct(k).match(tracker,1:2)
            s(obj_count,i)={[aux_struct(i).match(j,1:2);aux_struct(k).match(tracker,1:2)]};
            prev=aux_struct(k).match(tracker,1:2);
        else
            obj_count=obj_count+1;
        end
        
    end
end

