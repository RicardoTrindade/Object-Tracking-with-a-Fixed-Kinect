function [rgb,cmap,c] = RGBFormsMatching(xyz,imrgb,R_d_to_rgb,T_d_to_rgb,RGB_cam,prev_c,psize)

%Getting the RGB image with the heads only
rgb = get_rgbd(xyz,imrgb,R_d_to_rgb,T_d_to_rgb,RGB_cam);

%Obtaining the centroids in the rgb image
aux = rgb2gray(rgb);
pos = aux~=0;
aux(pos) = 1;
BW_aux = bwlabel(aux);
BW_aux = medfilt2(BW_aux,[5,5]);
s = regionprops(BW_aux,'centroid');
c = cat(1, s.Centroid);
c_filt = find(isnan(c));
[l,~] = ind2sub(size(c),c_filt);
c(l,:) = [];

%Normalizing the centroids coordinates
size_aux = size(BW_aux,2);
prev_c = prev_c.*(size_aux/psize);

%Cost Matrix
N = size(prev_c,1);
M = size(c,1);
m_cost = zeros(M,N);
for i=1:N
    for j=1:M
        m_cost(j,i) = prev_c(i,1)-c(j,2);
    end
end
[cmap,~,~]=assignDetectionsToTracks(m_cost,1000,1000);

end