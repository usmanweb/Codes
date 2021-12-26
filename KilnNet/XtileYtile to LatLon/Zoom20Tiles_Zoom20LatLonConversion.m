clc;
clear;
close all;
pkg load image
path='attique/';
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
end


lat=((Xtile-368247)*0.000686646)+72.85549164;
lon=((Ytile-368247)*0.000686646)+72.85549164;

T = [lon lat];

for i=1:length(names)
  images{i} = imread([path names{i}]);
  filename{i}=strcat(num2str(lon(i),16),'_',num2str(lat(i),16),'.jpg');
  imwrite(images{i},filename{i});
end

%filename = 'KilnsLatLon.xls'
%writetable(T,filename,'Sheet',1);
%xlswrite(filename, T)

%images = cell(length(names),1);
%% To remove empty cells
% images = images(~cellfun('isempty',images));
%for i=1:length(names) #15000
    %disp(i)
    %images{i} = imread([path names{i}]);
    %slidingwindow(images{i},StartCornerLat(i),StartCornerLon(i),diff);
%end


