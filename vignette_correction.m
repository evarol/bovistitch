function dataVolume=vignette_correction(dataVolume,numbins,numiter,templatetype,visual,slice)

vec=@(x)(x(:));
if nargin<2
    numbins=20;
    numiter=1;
    templatetype='middle_slice';
    visual=0;
end
if nargin<3
    numiter=1;
    templatetype='middle_slice';
    visual=0;
end
if nargin<4
    templatetype='middle_slice';
    visual=0;
end
if nargin<5
    visual=0;
end

if visual==1
    dataVolume0=dataVolume;
end
for k=1:numiter
    Xh=reshape(dataVolume,[size(dataVolume,1)],[]);
    
    Xht=zeros(size(Xh));
    if strcmpi(templatetype,'middle_slice');
        template=Xh(1024,:);
    elseif strcmpi(templatetype,'middle_20_slices');
        template=vec(Xh(1024-10:1024+10,:))';
    elseif strcmpi(templatetype,'random');
        r=randperm(numel(Xh));
        template=Xh(r(1:size(Xh,2)));
    end
    tic
    for i=1:size(Xh,1)
        Xht(i,:)=linhistmatch(Xh(i,:),template,numbins,'regular');
        clc
        fprintf(['Slice ' num2str(slice) ' Horizontal vignette correction progress (round ' num2str(k) '/' num2str(numiter) '):\n']);
        fprintf(['\n' repmat('.',1,50) '\n\n'])
        for tt=1:round(i*50/size(Xh,1))
            fprintf('\b|\n');
        end
        T=toc;
        disp(['Time elapsed (minutes): ' num2str(T/60) ' Time remaining (minutes): ' num2str((size(Xh,1)-i)*(T/i)*(1/60))]);
    end
    
    dataVolume=reshape(Xht,size(dataVolume));
    
    if visual==1
        imagesc([imgaussfilt(max(dataVolume0,[],3),200) imgaussfilt(max(dataVolume,[],3),200)]);title('Left: uncorrected vignetting, Right: corrected vignetting');colormap(gray(256));drawnow
    end
    Xv=reshape(permute(dataVolume,[2 1 3]),[2048],[])';
    Xvt=zeros(size(Xv));
    if strcmpi(templatetype,'middle_slice');
        template=Xv(:,1024);
    elseif strcmpi(templatetype,'middle_20_slices');
        template=vec(Xv(:,1024-10:1024+10));
    elseif strcmpi(templatetype,'random');
        r=randperm(numel(Xv));
        template=Xv(r(1:size(Xv,2)));
    end
    tic
    for i=1:size(Xv,2)
        Xvt(:,i)=linhistmatch(Xv(:,i),template,numbins,'regular');
        clc
        fprintf(['Slice ' num2str(slice) ' Vertical vignette correction progress (round ' num2str(k) '/' num2str(numiter) '):\n']);
        fprintf(['\n' repmat('.',1,50) '\n\n'])
        for tt=1:round(i*50/size(Xv,2))
            fprintf('\b|\n');
        end
        T=toc;
        disp(['Time elapsed (minutes): ' num2str(T/60) ' Time remaining (minutes): ' num2str((size(Xv,2)-i)*(T/i)*(1/60))]);
    end
    
    tmp=reshape(Xvt',size(dataVolume));
    
    dataVolume=permute(tmp,[2 1 3]);
    
    if visual==1
        imagesc([imgaussfilt(max(dataVolume0,[],3),200) imgaussfilt(max(dataVolume,[],3),200)]);title('Left: uncorrected vignetting, Right: corrected vignetting');colormap(gray(256));drawnow
    end
end


end


function [atransform,distance]=linhistmatch(a,b,nbins,type)
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
abins=quantile(a,nbins)';
bbins=quantile(b,nbins)';
ahat=histc(a,abins)';
bhat=histc(b,bbins)';
% weighted linear regression of the matching quantilesy
if strcmpi(type,'non-negative');
    beta=lsqnonneg(diag(bhat)*[abins ones(size(abins,1),1)],diag(bhat)*bbins);
elseif strcmpi(type,'regular');
    beta=linsolve(diag(bhat)*[abins ones(size(abins,1),1)],diag(bhat)*bbins);
end
%transformed time series with nan's put back in
atransform=nan(size(a_nan_idx));
atransform(a_nan_idx)=a*beta(1) + beta(2);

%wasserstein distance computation
% atransformhat = histc(atransform,bbins)'; %discretizing the transformed time series
% distance=wdist(atransformhat,bhat,1); %wasserstein distance computation between atransform and b
distance=[];
end
