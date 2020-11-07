clear all
clc
close all

%% parameters
numbins=200; %% a "pseudo" parameter -- defines the fineness of quantile estimation
numiter=5; %% it converges eventually
template='middle_slice'; %which portion of the FOV to normalize to
visual=1; %to see updates
downsample_factor=10; %% just to speed up vignette estimation (reduces X and Y resolutions)
downslice_factor=20; %% just to speed up vignette estimation (use only a random subset of slices)
channel=1; %which channel to stitch/correct
%% adding dependencies
disp('Select imarisReader folder...');
selpath = uigetdir('[]','Select imarisReader folder');
addpath(genpath(selpath));

%% import files
disp('Select ims files to stitch...');
[file,path] = uigetfile('*.ims','Select ims files to stitch','MultiSelect','on');


%% vignette correction - downsampling for speed up (don't need every single pixel to estimate the vignette field)

fprintf(['Loading slices for vignette estimation...\n']);
tic

t=0;
for i=1:length(file)
    imsObj=ImarisReader([path file{i}]);
    if i==1
        r=randperm(imsObj.DataSet.SizeZ)-1;
        slicesamples=r(1:ceil(length(r)/downslice_factor));
        dataVolume=zeros([round(imsObj.DataSet.SizeX/downsample_factor) round(imsObj.DataSet.SizeY/downsample_factor) length(file)*length(slicesamples)]);
    end
    for s=1:length(slicesamples)
        t=t+1;
        dataVolume(:,:,t)=imresize(imsObj.DataSet.GetDataSlice(slicesamples(s),channel,0),[round(imsObj.DataSet.SizeX/downsample_factor) round(imsObj.DataSet.SizeY/downsample_factor)]);
        
        clc
        fprintf(['Loading slices in dataset (' num2str(i) '/' num2str(length(file)) ')...\n']);
        fprintf(['\n' repmat('.',1,50) '\n\n'])
        for tt=1:round(s*50/length(slicesamples))
            fprintf('\b|\n');
        end
        T=toc;
        disp(['Time elapsed (minutes): ' num2str(T/60) ' Time remaining (minutes): ' num2str((length(slicesamples)-s)*(T/s)*(1/60))]);
    end
end
fprintf(['Loading slices...done\n']);

dataVolume=double(dataVolume);

%%%%%%% main solver %%%%%%%%
[~,vfield,vfield_corrected,S,D]=vignette_correction(dataVolume,numbins,numiter,template,visual,0); %% the main meat of the code - vfield is the smoothed vignette field, S and D are the multiplicative / additive corrections, respectively.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% upsampling (want to vignette correct the full res images)
for k=1:size(S,1)
    S{k,1}=imresize(S{k,1},[imsObj.DataSet.SizeX 1]);
    S{k,2}=imresize(S{k,2},[1 imsObj.DataSet.SizeY]);
    D{k,1}=imresize(D{k,1},[imsObj.DataSet.SizeX 1]);
    D{k,2}=imresize(D{k,2},[1 imsObj.DataSet.SizeY]);
end
clear dataVolume

%% stitching



tic
for s=1:imsObj.DataSet.SizeZ
    slice=s-1;
    
    disp('Loading slices...');
    for i=1:length(file)
        imsObj=ImarisReader([path file{i}]);
        
        dataVolume(:,:,i)=double(imsObj.DataSet.GetDataSlice(slice,channel,0));
        dataVolume_corrected(:,:,i)=double(dataVolume(:,:,i));
        for k=1:size(S,1)
            dataVolume_corrected(:,:,i)=dataVolume_corrected(:,:,i).*S{k,1} + D{k,1}; %% Vignette correction applied
            dataVolume_corrected(:,:,i)=dataVolume_corrected(:,:,i).*S{k,2} + D{k,2};
        end
        position(:,i)=[imsObj.DataSet.ExtendMinX imsObj.DataSet.ExtendMaxX imsObj.DataSet.ExtendMinY imsObj.DataSet.ExtendMaxY];
    end
    disp('Loading slices...(done)');
    
    
    
    
    
    disp(['Stitching slice ' num2str(slice) ]);
    position(1:2,:)=position(1:2,:)-min(position(1,:));
    position(3:4,:)=position(3:4,:)-min(position(3,:));
    position=round(position+1);
    counts=zeros(max(position(2,:)),max(position(4,:)));
    I=zeros(max(position(2,:)),max(position(4,:)));
    Icorr=zeros(size(I));
    vignette_field=zeros(size(I));
    vignette_field_corrected=zeros(size(I));
    for i=1:size(dataVolume,3)
        I(position(1,i):position(2,i),position(3,i):position(4,i))=I(position(1,i):position(2,i),position(3,i):position(4,i))+imresize(dataVolume(:,:,i),size(I(position(1,i):position(2,i),position(3,i):position(4,i))));
        Icorr(position(1,i):position(2,i),position(3,i):position(4,i))=Icorr(position(1,i):position(2,i),position(3,i):position(4,i))+imresize(dataVolume_corrected(:,:,i),size(I(position(1,i):position(2,i),position(3,i):position(4,i))));
        vignette_field(position(1,i):position(2,i),position(3,i):position(4,i))=vignette_field(position(1,i):position(2,i),position(3,i):position(4,i))+imresize(vfield,size(I(position(1,i):position(2,i),position(3,i):position(4,i))));
        vignette_field_corrected(position(1,i):position(2,i),position(3,i):position(4,i))=vignette_field_corrected(position(1,i):position(2,i),position(3,i):position(4,i))+imresize(vfield_corrected,size(I(position(1,i):position(2,i),position(3,i):position(4,i))));
        counts(position(1,i):position(2,i),position(3,i):position(4,i))=counts(position(1,i):position(2,i),position(3,i):position(4,i))+1;
    end
    
    I=I./counts; %% Raw stitch output
    Icorr=Icorr./counts; %% Bovi-stitch output
    
    vignette_field=vignette_field./counts; %% Raw vignette field
    vignette_field_corrected=vignette_field_corrected./counts; %% Raw vignette field
    
    volume_raw(:,:,slice+1)=I; %% output stitched volume
    volume_corrected(:,:,slice+1)=Icorr; %% output bovi-stitched volume
    
            clc
        fprintf(['Stitching slice ' num2str(slice) '\n']);
        fprintf(['\n' repmat('.',1,50) '\n\n'])
        for tt=1:round(s*50/imsObj.DataSet.SizeZ)
            fprintf('\b|\n');
        end
        T=toc;
        disp(['Time elapsed (minutes): ' num2str(T/60) ' Time remaining (minutes): ' num2str((imsObj.DataSet.SizeZ-s)*(T/s)*(1/60))]);
    
    disp(['Stitching slice ' num2str(slice) ' (done)']);
end

%% save results
disp('All slices finished');
disp('Saving stitched + vignetted corrected volume');
save([path 'stitched_raw_volume.mat'],'volume_raw','vignette_field','vfield');
save([path 'bovistitched_volume.mat'],'volume_corrected','vignette_field_corrected','vfield_corrected');
disp('Saving stitched + vignetted corrected volume (done)');