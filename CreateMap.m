function [ s ] = CreateMap( tracks )


n=length(tracks);
s=cell(n-1,1);
for i=1:n-1
    u=tracks(i).match;
    v=tracks(i+1).match;
    if isempty(v)
        continue;
    end
    if isempty(u)
        continue;
    end
    m=max(max(u(:)),max(v(:)));
    aux=zeros(m,1);
    for j=1:size(v(:,1),1)
      z=v(j,1);
      aux(z)=v(j,2);
    end
    s(i)={aux'};
end
end