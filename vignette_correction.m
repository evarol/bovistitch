function [dataVolume,vfield,vfield_corrected,S,D]=vignette_correction(dataVolume,numbins,numiter,templatetype,visual,slice)

vec=@(x)(x(:));
if nargin<2
    slice=0;
    numbins=20;
    numiter=1;
    templatetype='middle_slice';
    visual=0;
end
if nargin<3
    slice=0;
    numiter=1;
    templatetype='middle_slice';
    visual=0;
end
if nargin<4
    slice=0;
    templatetype='middle_slice';
    visual=0;
end
if nargin<5
    slice=0;
    visual=0;
end

% if visual==1
    dataVolume0=dataVolume;
% end

for k=1:numiter
    Xh=reshape(dataVolume,[size(dataVolume,1)],[]);
    
    Xht=zeros(size(Xh));
    if strcmpi(templatetype,'middle_slice');
        template=Xh(round(size(dataVolume,1)/2),:);
    elseif strcmpi(templatetype,'middle_20_slices');
        template=vec(Xh(round(size(dataVolume,1)/2)-10:round(size(dataVolume,1)/2)+10,:))';
    elseif strcmpi(templatetype,'random');
        r=randperm(numel(Xh));
        template=Xh(r(1:size(Xh,2)));
    end
    tic
    for i=1:size(Xh,1)
        [Xht(i,:),~,Bh(i,:)]=linhistmatch(Xh(i,:),template,numbins,'regular');
        clc
        fprintf(['Slice ' num2str(slice) ' Horizontal vignette correction progress (round ' num2str(k) '/' num2str(numiter) '):\n']);
        fprintf(['\n' repmat('.',1,50) '\n\n'])
        for tt=1:round(i*50/size(Xh,1))
            fprintf('\b|\n');
        end
        T=toc;
        disp(['Time elapsed (minutes): ' num2str(T/60) ' Time remaining (minutes): ' num2str((size(Xh,1)-i)*(T/i)*(1/60))]);
    end
    
    S{k,1}=Bh(:,1);D{k,1}=Bh(:,2);
    
    dataVolume=reshape(Xht,size(dataVolume));
    
    if visual==1
        imagesc([imgaussfilt(max(dataVolume0,[],3),size(dataVolume,1)/5) imgaussfilt(max(dataVolume,[],3),size(dataVolume,1)/5)]);title('Left: uncorrected vignetting, Right: corrected vignetting');colormap(gray(256));drawnow
    end
    Xv=reshape(permute(dataVolume,[2 1 3]),[size(dataVolume,2)],[])';
    Xvt=zeros(size(Xv));
    if strcmpi(templatetype,'middle_slice');
        template=Xv(:,round(size(dataVolume,2)/2));
    elseif strcmpi(templatetype,'middle_20_slices');
        template=vec(Xv(:,round(size(dataVolume,2)/2)-10:round(size(dataVolume,2)/2)+10));
    elseif strcmpi(templatetype,'random');
        r=randperm(numel(Xv));
        template=Xv(r(1:size(Xv,2)));
    end
    tic
    for i=1:size(Xv,2)
        [Xvt(:,i),~,Bv(:,i)]=linhistmatch(Xv(:,i),template,numbins,'regular');
        clc
        fprintf(['Slice ' num2str(slice) ' Vertical vignette correction progress (round ' num2str(k) '/' num2str(numiter) '):\n']);
        fprintf(['\n' repmat('.',1,50) '\n\n'])
        for tt=1:round(i*50/size(Xv,2))
            fprintf('\b|\n');
        end
        T=toc;
        disp(['Time elapsed (minutes): ' num2str(T/60) ' Time remaining (minutes): ' num2str((size(Xv,2)-i)*(T/i)*(1/60))]);
    end
    
    S{k,2}=Bv(1,:);D{k,2}=Bv(2,:);
    tmp=reshape(Xvt',size(dataVolume));
    
    dataVolume=permute(tmp,[2 1 3]);
    
    if visual==1
        imagesc([imgaussfilt(max(dataVolume0,[],3),size(dataVolume,1)/5) imgaussfilt(max(dataVolume,[],3),size(dataVolume,1)/5)]);title('Left: uncorrected vignetting, Right: corrected vignetting');colormap(gray(256));drawnow
    end
end

vfield=imgaussfilt(max(dataVolume0,[],3),size(dataVolume,1)/5);
vfield_corrected=imgaussfilt(max(dataVolume,[],3),size(dataVolume,1)/5);
end


function [atransform,distance,beta]=linhistmatch(a,b,nbins,type)
%takes as input two traces a (moving), b (reference) and outputs normalized
%trace atransform that has a similar histogram as b
%input : a - moving time series
%        b - reference time series
%        nbins - number of bins to discretize both time series into
%output: atransform - moved time series
%        distance - histogram distance between atransform and b

% discarding nans from time series
a_nan_idx=~isnan(a);
b_nan_idx=~isnan(b);
a=a(a_nan_idx);
b=b(b_nan_idx);

% discretizing time series using quantiles
abins=quantile(a,linspace(0,1,nbins))'; %<---FOR BOVEY: changed histograms to include the max and min, before quantile(a,nbins) DID NOT include min and max!
bbins=quantile(b,linspace(0,1,nbins))';
% weighted linear regression of the matching quantiles
if strcmpi(type,'non-negative')
    beta=lsqnonneg([abins ones(size(abins,1),1)],bbins); %<---FOR BOVEY: also changed the regression slightly, the weights are not used anymore since all bins are equal weights (if we use a non-linear bins, then weights may be used)
elseif strcmpi(type,'regular')
    beta=linsolve([abins ones(size(abins,1),1)],bbins);
end
%transformed time series with nan's put back in

atransform=nan(size(a_nan_idx));
atransform(a_nan_idx)=a*beta(1) + beta(2);

%wasserstein distance computation
% atransformhat = histc(atransform,bbins)'; %discretizing the transformed time series
% distance=wdist(atransformhat,bhat,1); %wasserstein distance computation between atransform and b
distance=[]; %% turned this off since it seems to be not too useful
end
