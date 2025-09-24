function [bb,bbabs]=pinv_recon(data,wfn,varargin)
%% Reconstruction with pseudoinverse of encoding matrix
% Reconstructing images with Pinv-Recon. Looks for waveform name saved as
% [wfn_F(mtx_reco(1))(fnstr)(_b0)(_sens).h5] if the encoding matrix size 
% > 100MB.
%
% Inputs      data  raw data in (ncoils,npts,ntimesteps,nslices,nspec)             
%              wfn  waveform name (struct also accepted - requires wf.k and
%                   wf.t and will not read or save reconstruction matrices)       
%
% Optional inputs
%            shift  spatial shift in x, y ,z                      ([0 0 0])
%      debug_slice  choose only a single slice to recon (debugging)    ([])
%         mtx_reco  reconstructed image matrix size                    ([])
%               b0  B0 map for off-res correction (x,y,z)              ([])
%             sens  Sensitivity map for coil sense encoding            ([])
%                   (x,y,z,ncoils)
%             mode  decomposition method ('svd','qr','eig',    ('cholesky')
%                   'cholesky')
%           thresh  threshold for regularization           (default depends 
%                                                                  on mode)
%     saveSRFNoise  saves SRF and Noise maps where data is loaded       (0)
%            fnstr  custom string in the saving file name              ([])
% Outputs
%               bb  complex image
%            bbabs  absolute image (root sum of squares coil combine)
%
% 08/2024   Kylie Yeung

if (nargin<1), help(mfilename); return; end
tStart=tic;

%% reading waveform file
if isstruct(wfn)
    wf = wfn;
    wfn='unknown';Fdir='unknown';name='unknown';
    warning('Waveform provided in struct format. Will not be attempting to read or write reconstruction matrix.')
else
        % load waveform structure from fidall.e
    [Fdir,name,ext]=fileparts(wfn);
    wfn=[Fdir '/' name '.mat'];
    if ~exist(wfn,'file'), error('file not found: wfn=%s',wfn); end
    try
        wf = load_waveform(wfn,true,true);       % load waveform
    catch
        load(wfn)
    end
end

%% parsing inputs
[mtx_acq,mtx_reco,shift,b0,sens,debug_slice,mode,thresh,useGPU,saveSRFNoise,fnstr]=parseInputs(varargin,wf);

%% getting filename suffix
fnstr=[fnstr num2str(mtx_reco(1))];
if ~isempty(b0);fnstr=[fnstr '_b0'];end
if ~isempty(sens);fnstr=[fnstr '_sens'];end

if saveSRFNoise && exist([Fdir '/SRFNoise' fnstr '.mat'],'file');saveSRFNoise=0;disp('SRF and noise maps already exist.');end

%% ensure mtx_reco is a 1x3
[ncoils,npts,ntimesteps,nslices,nspec]=size(data);
if length(mtx_acq)==1; mtx_acq=[mtx_acq mtx_acq nslices];end
if length(mtx_acq)==2; mtx_acq=[mtx_acq nslices];end
if length(mtx_reco)==1; mtx_reco=[mtx_reco mtx_reco nslices];end
if length(mtx_reco)==2; mtx_reco=[mtx_reco nslices];end

%% get reconstruction matrix
k=wf.k;
if ~isreal(k)
    if size(k,1)~=1, error('size(k,1)~=1'); end
    k_temp=k;
    k(1,:,1) = real(k_temp);
    k(1,:,2) = imag(k_temp);
    clear k_temp
end

% modify to dim=2 if the problem is separable
if isfield(wf,'kz')
    separable3d = true;
    [F_z,k,dim] = getFz(k,mtx_acq, mtx_reco,npts);
    wf.t=wf.t(:,1:npts/mtx_acq(3));
else
    separable3d=false; dim=size(k,3);
end

% determine encoding matrix size
if isempty(sens)
    memE = size(k,2)*prod(mtx_reco(1:2))*4*2*1d-6;
    fprintf('size(E) = [%d %d] = %g [MB]\n',size(k,2),prod(mtx_reco(1:2)),memE);
else
    memE = size(k,2)*ncoils*prod(mtx_reco(1:2))*4*2*1d-6;
    fprintf('size(E) = [%d %d] = %g [MB]\n',size(k,2)*ncoils,prod(mtx_reco(1:2)),memE);
end
% Determine the filename without extensions
Fn = [Fdir '/' name '_F' fnstr '.h5'];

readfile=0;savefile=0;
if ~exist(Fn, 'file'); readfile=0; savefile=1; disp('File not found. Not reading.');end
if memE < 100; readfile=0; savefile=0; disp('Small matrix. Not reading or saving.');end

% Check if the reconstruction matrix should be calculated
if readfile
    fprintf('Loading recon matrix file=\n\t%s\n', Fn);
    loadStart = tic;
    F = readh5complex(Fn); % Load waveform
    fprintf('Reading recon matrix: %g s\n', toc(loadStart));
end

% Check if the encoding matrix should be calculated
if ~readfile || saveSRFNoise
    E = getE(k, mtx_acq, mtx_reco, dim, shift); % get gradient only encoding matrix

    % modify encoding matrix with B0 encoding
    if ~isempty(b0)
        switch ndims(b0)
            case 2
                b0=imresize(b0,mtx_reco(1:2));
            case 3
                b0=imresize3(b0,mtx_reco);
        end

        for i=1:length(shift); b0=circshift(b0,-shift(i),i);end
        if ~isempty(debug_slice);b0=b0(:,:,debug_slice);end
        [E] = getE_b0(E,b0,wf.t);
    end
    
    % modify encoding matrix with coil sensitivity encoding
    if ~isempty(sens)
        sens=imresizen(sens,[mtx_reco./size(sens,1,2,3) 1]);
        for i=1:length(shift); sens=circshift(sens,-shift(i),i);end
        if ~isempty(debug_slice);sens=sens(:,:,debug_slice,:);end
        [E] = getE_sens(E,sens);
    end
end
if ~isempty(debug_slice);mtx_reco(3)=1;end
if ~readfile
    % obtain reconstruction matrix
    if ~isempty(debug_slice);mtx_reco(3)=1;end
    if ~iscell(E)
       F = getF(E, thresh, mode, useGPU);
    else
       [F] = getF_slice(E,thresh,mode,useGPU);
    end

    % Save reconstruction matrix if memory threshold allows
    if savefile
        savestart = tic;
        saveh5complex(Fn, gather(F));
        fprintf('Saving recon matrix: %g s\n', toc(savestart));
    end
end

%% obtain SRF and noise maps if required
if saveSRFNoise
   if iscell(E);E=cat(3,E{:});end
   SRF = getSRF(E, F, dim, mtx_reco);
   NoiseMat = getNoiseMat(F, [], dim);
   save([Fdir '/SRFNoise' fnstr '.mat'],"SRF","NoiseMat")

   clear SRF NoiseMat
end

clear E saveSRFNoise

%% reconstruction
if separable3d
    fprintf('Pinv in z... \n')
    nslices = size(F_z,1); nz_acq = size(F_z,2);
    data = reshape(permute(data,[1 3 4 5 2]),[ncoils*ntimesteps*size(F,2),nz_acq]);
    data = F_z*data.';
    data = permute(reshape(data,[nslices,ncoils,ntimesteps,size(F,2)]),[2 4 3 1 5]);
    if ~isempty(debug_slice); data=data(:,:,:,debug_slice);end
end

if length(mtx_reco)==2;mtx_reco(3)=nslices;end
if ~isempty(sens);ncoils=1;end
if ndims(F)==3
    % slicewise reconstruction, required for B0 and sense recon
    bb=zeros(mtx_reco(1),mtx_reco(2),mtx_reco(3),ntimesteps,ncoils);
    for ls=1:mtx_reco(3) %round(mtx_reco(3)/2) %
        dd_s=reshape(permute(data(:,:,:,ls),[2 1 3]),[size(F,2),ntimesteps*ncoils]);
        data_inverted=F(:,:,ls)*dd_s;
        bb(:,:,ls,:,:)=reshape(data_inverted,[mtx_reco(1),mtx_reco(2),1,ntimesteps,ncoils]);
    end
else
    dd=reshape(permute(data,[2 4 3 1 5]),[size(F,2),mtx_reco(3)*ntimesteps*ncoils*nspec]);
    data_inverted=F*dd;
    bb=reshape(data_inverted,[mtx_reco(1),mtx_reco(2),mtx_reco(3),ntimesteps,ncoils,nspec]);
end
clear F data

%% coil combine
if ncoils>1
    bbabs = sqrt(mean(bb.*conj(bb),5));
else
    bbabs = abs(bb);
end

%% timing
times.total_time=toc(tStart);
try 
    fprintf('recon_pinv: runtime = %.2f [s] (%.2f [s] for svd) \n ',times.total_time, times.svd_time(1))
catch
    fprintf('recon_pinv: runtime = %.2f [s] \n ',times.total_time)
end

end

function [mtx_acq,mtx_reco,shift,b0,sens,debug_slice,mode,thresh,useGPU,saveSRFNoise,fnstr]=parseInputs(varargin,wf)
if ~isfield(wf,'mtx') && ~isfield(wf,'npix');error('Acquisition matrix size is required in wf.mtx or wf.npix');end
try mtx_acq=wf.mtx; catch; mtx_acq=wf.npix; end
if length(mtx_acq)==1;mtx_acq=[mtx_acq mtx_acq];end
mtx_reco=mtx_acq;
shift=[0 0 0];
b0=[];
sens=[];
debug_slice=[];
mode='cholesky';
thresh=[];
useGPU=canUseGPU;
saveSRFNoise=0;
fnstr=[];

% Parse varargin
for i = 1:2:length(varargin)
    switch lower(varargin{i})
        case 'mtx_reco'% make a filename suffix for checking if the ones with different options
% are chosen
            mtx_reco = varargin{i+1};
        case 'shift'
            shift = varargin{i+1};
        case 'b0'
            b0 = varargin{i+1};
        case 'sens'
            sens = varargin{i+1};
        case 'debug_slice'
            debug_slice = varargin{i+1};
        case 'mode'
            mode = varargin{i+1};
        case 'thresh'
            thresh = varargin{i+1};
        case 'usegpu'
            useGPU = varargin{i+1};
        case 'savesrfnoise'
            saveSRFNoise = varargin{i+1};
        case 'fnstr'
            fnstr = varargin{i+1};
        otherwise
            error('Unknown parameter %s', varargin{i});
    end
end

end