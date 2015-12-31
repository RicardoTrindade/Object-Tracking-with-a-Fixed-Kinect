function [ output,v ] = BuildTrajectory( tracks )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
output=struct();
n=length(tracks);

for i=1:n
    while(~isempty(tracks(i).match))
       aux=tracks(i).match(i,2);
       v=aux;
       tracks(i).match(i,:)=[];
       for k=i+1:n
           if(tracks(k).match(k,2)==aux)
              v=[v;tracks(k).match(k,2)];      
           end
       end
       
    end
end





end

