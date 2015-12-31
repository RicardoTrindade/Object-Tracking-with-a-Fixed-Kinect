function [ tracked_objs ] = peopletracking( file_names, Depth_cam, RG_cam, Rdtrgb,Tdtrgb )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
nfiles=length(file_names);
loadvec=round(nfiles*(rand(1,30)));
cut=all(loadvec,1);
loadvec=loadvec(:,cut~=0);
depth = zeros(480,640,length(loadvec));

%This code chunk takes random images from the whole dataset and
%computes the median of them all for using in background removal
for i=1:length(loadvec)
    load(file_names(i).depth)
    depth(:,:,i) = depth_array;
end

med = median(depth,3); %The Brackground Median that will be subtracted to all images
medreshaped = reshape(med,480*640,1);
ground = get_xyzasus(medreshaped,[480 640],find(medreshaped(:)>0),Depth_cam.K,1,0);%Points that correspond to the background



%Obtaing the normal vector to the ground, the position of the ground and its points through RANSAC
[normal_vec,d,no_groundpoints,ground_points] = RANSAC(ground,3,500,0.02,0.2);
% Orientation of the normal vector to the ground
if(normal_vec(1)<0)
    normal_vec(1) = -normal_vec(1);
end
if(normal_vec(2)>0)
    normal_vec(2) = -normal_vec(2);
end
if(normal_vec(3)>0)
    normal_vec(3) = -normal_vec(3);
end
if(d>0)
    d = -d;
end

nx = cross(normal_vec,[0 0 1]);%Obtaining the rotation vector
nx = nx/norm(nx);
theta = dot(normal_vec,[0 0 1])/(norm([0 0 1])*norm(normal_vec));
theta = acos(theta);

%Obtaining the Rotation Matrix using Rodriguez's Formula
X = [0 -nx(3) nx(2);nx(3) 0 -nx(1);-nx(2) nx(1) 0];
R = eye(3)+sin(theta)*X+(1-cos(theta))*X^2;

%Rotate the Ground
rotate = R*ground';
range_x = round([min(rotate(1,:))-1 1+max(rotate(1,:))]); %Used for building the image grill later
range_y = round([min(rotate(2,:))-1 1+max(rotate(2,:))]);
clear ground X theta nx ground_points x_mean y_mean z_mean filename2
clear depth depth_array filename frame_date frame_toc rgb_array loadvec i
clear x_m y_m z_m
%%
load(file_names(1).depth); %Load the image 1
[present_matrix,area1,prev_BW,prev_num] = first_procedure(depth_array,med,Depth_cam,R,d,range_x,range_y);
tracks = struct();
Forms=cell(1);
%2) Track the objects in the remaining frames
figure(1);
for z=2:nfiles
    load(file_names(z).depth); %Load the image z
    imrgb = imread(file_names(z).rgb); %Load the rgb image
    
    %Removal of the background on the image
    im2 = double((depth_array));
    mask2 = med-im2>800;
    im_noback2 = reshape(mask2.*im2,480*640,1);
    
    %Rotation of the image with no background
    xyz2 = get_xyzasus(im_noback2,[480 640],find(im_noback2(:)>0),Depth_cam.K,1,0);
    rotate2 = R*xyz2';
    aux = rotate2;
    
    %Select only the points that are above 1.4m in order to select the heads
    altura2 = find(rotate2(3,:)<-0.8+d);
    aux(:,altura2) = 0;
    im_new = aux'*R;
    rotate2(:,altura2) = [];
    
    %Creating a 2D grid of the 3D points
    delta = 0.05; %Grid division
    u2 = zeros((range_x(2)-range_x(1))/delta,(range_y(2)-range_y(1))/delta);%Grid
    
    %Selection of xy coords. only
    x2 = rotate2(1:2,:);
    cut = all(x2,1);
    s2 = x2(:,cut~=0);
    
    %Division of coords by grid division
    ptsfix2 = fix(s2/delta);
    clear cut rotate2 rotate xyz2
    
    %Representation of the points
    for j=1:size(ptsfix2,2)
        u2(ptsfix2(1,j)+abs(range_x(1))/delta+1,ptsfix2(2,j)+abs(range_y(1))/delta+1) = u2(ptsfix2(1,j)+abs(range_x(1))/delta+1,ptsfix2(2,j)+abs(range_y(1))/delta+1)+1;
    end
    
    %Remotion of noise from the grid
    u2 = medfilt2(u2,[5 5]);
    [BW2,num2] = bwlabel(u2,8);
    [BW2,num2] = AreaFilt(BW2,num2);
%     imagesc(BW2);
    Forms{end+1}=BW2;
    %Joint Blobs Segmentation (to separate different forms)
    % blobs_only=nonzeros(BW2);
    % bigger_blob=mode(blobs_only);
    %
    % if sum(blobs_only==bigger_blob)>65
    %     D=bwdist(~BW2);
    %     D = -D;
    %     D(~BW2) = -Inf;
    %     L=watershed(D);
    %     BW2(L == 0) = 0;
    %     [BW2,num2]=bwlabel(BW2,8);
    % end
    
    %Obtaining the Cost Matrix for the Hungarian Method
    area2 = formsArea(BW2,num2); %Determine the area of all the forms in the frame
    s = regionprops(BW2,'centroid'); %Obtain the centroids of all the forms in the frame
    future_matrix = cat(1, s.Centroid); %Refresh future_matrix
    %Checking if the frame is empty
    if isempty(future_matrix)
        area2 = area1;
        future_matrix = present_matrix;
        continue
    end
    c = CostMatrixAreas(area2,area1,num2);%Cost Matrix for the Areas
    cs = CostMatrixShape(BW2,num2,prev_BW,prev_num); %Cost Matrix for the Shape
    est_dist = zeros(size(present_matrix,1),size(future_matrix,1)); %Cost Matrix for the Distance of the centroids
    for i = 1:size(present_matrix, 1)
        diff = future_matrix - repmat(present_matrix(i,:),[size(future_matrix,1),1]);
        est_dist(i, :) = sqrt(sum(diff .^ 2,2));
    end
    est_dist = est_dist/16; %Divinding by the mean value of the max values
    est_dist(isnan(est_dist)) = 0;
    cost_matrix = est_dist+c+cs; %Sum of all cost Matrices
    [matches,unassignedtracks,unassigneddetections]=assignDetectionsToTracks(cost_matrix,2.5,2.5); %Hungarian Method for Objects Assignment
    
    %Save all the matches, unassigned objects, unassigned detections and
    %centroids of the objects
    tracks(z).match = matches;
    tracks(z).unassigned = unassignedtracks;
    tracks(z).undetect = unassigneddetections;
    tracks(z).place = s;
    
    %Obaining the RGB image
    [O_rgb,~,~] = RGBFormsMatching(im_new,imrgb,Rdtrgb,Tdtrgb,RG_cam.K,future_matrix,140);
    
    %Representing the results
    subplot(1,2,1);
    imshow(BW2);
    title('Forms on XY (no Tracking)');
    subplot(1,2,2);
    imagesc(O_rgb);
    title(sprintf('RGB image, frame=%g, num=%f',z,num2));
    pause(0.04);
    
    %Refreshing the Variables
    present_matrix = future_matrix;
    area1 = area2;
    prev_BW = BW2;
    prev_num = num2;
    
end

%% Tracking
%Representing the results

tracked_objs=cell(1);
clear map;
% map = Mapping(tracks);
[tracking,kappas] = BeginTrack(tracks);
count = 1;
for i=1:length(tracking)
    if i>1
        prev_track = tracking(i-1);
    else
        prev_track = 0;
    end
    if tracking(i)==prev_track
        count = count+1;
        map_aux = Mapping(tracks,tracking(i),kappas(i),count);
        map{1,i} = map_aux;
        map{2,i} = length(map_aux);
    else
        count = 1;
        map_aux = Mapping(tracks,tracking(i),kappas(i),count);
        map{1,i} = map_aux;
        map{2,i} = length(map_aux);
    end
end

figure(2);
start = tracking(1);
for i=start:nfiles
    hold on
    imshow(Forms{i});
    for j=1:size(map,2)
        aux = map{1,j};
        n_aux = map{2,j};
        if i<=n_aux && aux(i)~=0
                C = tracks(i).place(aux(i)).Centroid;
                text(C(1,1),C(1,2),num2str(j),'color','r');
        end
    end
    for i1 = 1:size(map,2)
        aux = map{1,i1};
        l = map{2,i1};
        count = 1;
        cent_aux = [];
        for j=1:l
            if aux(j)~=0
            C = tracks(j).place(aux(j)).Centroid;
            cent_aux(count,:) = [C(1,1)*delta+abs(range_x(1))*delta C(1,2)*delta+abs(range_x(2))*delta 0 j];
            count = count+1;
            end
        end
        tracked_objs{i1} = cent_aux;
    end
    title(sprintf('Tracking of People, Frame=%g',i));
    pause(0.05);
end

end

