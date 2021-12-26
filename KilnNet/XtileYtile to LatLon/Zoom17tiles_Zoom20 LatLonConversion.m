clc;
clear;
close all;
pkg load image
%path='Lahorez17/newtask/17/46317/';
path = 'Nepal/WholeKathmandu/Kathmandu2020_z17/newtask/17/all/';
file = dir(fullfile(path));
NF = length(file);
for k = 1 : NF
    thisfile = fullfile(path, file(k).name); 
end
%% Find Tile size from the names
names = arrayfun(@(file) [file.name], file, 'UniformOutput',false);
names{1}=[];
names{2}=[];
names = names(~cellfun('isempty',names));
[numbers]=regexp(names,'\d*','Match');

Xtile=zeros(length(names),1);
Ytile=zeros(length(names),1);
for i=1:length(numbers)
    Xtile(i)=str2double(numbers{i,1}(1));
    Ytile(i)=str2double(numbers{i,1}(2));
    %if i<length(numbers)
    %XtileSmall(i)=abs(str2double(numbers{i,1}(1))-str2double(numbers{i+1,1}(1)))/8;
    %YtileSamll(i)=abs(str2double(numbers{i,1}(2))-str2double(numbers{i+1,1}(2)))/8;
    %end
end

lon=((Xtile-42295)*0.0054931640625)+52.33612060546875;  % CENTRAL Lat Lon of each image tile for Zoom 17 
lat=((Ytile-42295)*0.0054931640625)+52.33612060546875;

%%  Converting zoom 17 tiles'names to zoom 20 lat lon

diff=6.866455078125E-4;  % this is the difference of lat/long in adjacent zoom 20 tiles.
clat=lat+4*diff;   % this is the CORNER lat for zoom 20 image from zoom 17 images
clon=lon-4*diff;     % this is the CORNER long for zoom 20 image from zoom 17 images

%% Corener tile centre
StartCornerLat=clat-0.5*diff;
StartCornerLon=clon+0.5*diff;


images = cell(length(names),1);
%% To remove empty cells
% images = images(~cellfun('isempty',images));
for i=1:length(names) #15000
    disp(i)
    images{i} = imread([path names{i}]);
    slidingwindow(images{i},StartCornerLat(i),StartCornerLon(i),diff);
    #pixelsSlidingWindow(images{i},StartCornerLat(i),StartCornerLon(i),diff);

end

