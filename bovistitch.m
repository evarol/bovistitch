clear all
clc
close all

%% parameters
numbins=50;
numiter=1;
template='random';
visual=1;
%% adding dependencies
selpath = uigetdir('Select imarisReader folder');
addpath(genpath(selpath));

%% import files

[file,path] = uigetfile('*.ims','MultiSelect','on');

slice=-1;

while 1==1
    slice=slice+1;
    try
        disp('Loading slices...');
        for i=1:length(file)
            imsObj=ImarisReader([path file{i}]);
            dataVolume(:,:,i)=imsObj.DataSet.GetDataSlice(slice,0,0);
            position(:,i)=[imsObj.DataSet.ExtendMinX imsObj.DataSet.ExtendMaxX imsObj.DataSet.ExtendMinY imsObj.DataSet.ExtendMaxY];
        end
         disp('Loading slices...(done)');
        %% vignette correction
        
        
        dataVolume=double(dataVolume);
        
        dataVolume_corrected=vignette_correction(dataVolume,numbins,numiter,template,visual,slice);
        
        
        
        %% stitching
        
        disp(['Stitching slice ' num2str(slice) ]);
        position(1:2,:)=position(1:2,:)-min(position(1,:));
        position(3:4,:)=position(3:4,:)-min(position(3,:));
        position=round(position+1);
        counts=zeros(max(position(2,:)),max(position(4,:)));
        I=zeros(max(position(2,:)),max(position(4,:)));
        Icorr=zeros(size(I));
        for i=1:size(dataVolume,3)
            I(position(1,i):position(2,i),position(3,i):position(4,i))=I(position(1,i):position(2,i),position(3,i):position(4,i))+imresize(dataVolume(:,:,i),size(I(position(1,i):position(2,i),position(3,i):position(4,i))));
            Icorr(position(1,i):position(2,i),position(3,i):position(4,i))=Icorr(position(1,i):position(2,i),position(3,i):position(4,i))+imresize(dataVolume_corrected(:,:,i),size(I(position(1,i):position(2,i),position(3,i):position(4,i))));
            counts(position(1,i):position(2,i),position(3,i):position(4,i))=counts(position(1,i):position(2,i),position(3,i):position(4,i))+1;
        end
        
        I=I./counts; %% Raw stitch output
        Icorr=Icorr./counts; %% Bovi-stitch output
        
        volume_raw(:,:,slice+1)=I; %% output stitched volume
        volume_corrected(:,:,slice+1)=Icorr; %% output bovi-stitched volume
        
        disp(['Stitching slice ' num2str(slice) ' (done)']);
    catch
        disp('slices finished');
        disp('Saving stitched + vignetted corrected volume');
        save([path 'stitched_raw_volume.mat'],'volume_raw');
        save([path 'bovistitched_volume.mat'],'volume_corrected');
        disp('Saving stitched + vignetted corrected volume (done)');
        break
    end
    
end