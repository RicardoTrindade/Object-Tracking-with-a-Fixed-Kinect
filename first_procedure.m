function [ present_matrix,area,BW,num ] = first_procedure( depth_array,med,Depth_cam,R,d,range_x,range_y )
%First_Procedure Builds the first matrix with the located objects
%Built for ease of code reading
im=double(depth_array);
mask=med-im>500;
im_noback=reshape(mask.*im,480*640,1);
xyz=get_xyzasus(im_noback,[480 640],find(im_noback(:)>0),Depth_cam.K,1,0);
rotate=R*xyz';
altura= rotate(3,:)<d+1.3;
rotate(:,altura)=[];
delta=0.05;
u1=zeros((range_x(2)-range_x(1))/delta,(range_y(2)-range_y(1))/delta);
x1=rotate(1:2,:);
cut=all(x1,1);
s1=(x1(:,cut~=0));
ptsfix1=fix(s1/delta);
for j=1:length(ptsfix1)
    u1(ptsfix1(1,j)+4/delta+1,ptsfix1(2,j)+1/delta+1)= u1(ptsfix1(1,j)+4/delta+1,ptsfix1(2,j)+1/delta+1)+1; %Atribution of points to grid position
    %This is way to specific, tem de ficar genérico
end

[BW,num]=bwlabel(u1);
s=regionprops(BW,'centroid');
present_matrix=cat(1, s.Centroid);
area=formsArea(BW,num);


end

