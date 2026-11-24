function RMSE = main_measured_rmse(dataFile)
% Figure 4. Read the NWPU recording from this folder or from the supplied MAT file.
algorithmDir = fileparts(mfilename('fullpath'));
if nargin == 0
    dataFile = fullfile(algorithmDir, 'ACDataset_NWPU_Aug2026.mat');
end
if ~isfile(dataFile)
    error('DOA:MissingRecording', ...
        'Recording not found: %s. See README.md for the download link.',dataFile);
end
[~,fileInfo] = fileattrib(dataFile);
dataFile = fileInfo.Name;

snapshot = 10;                      % Snapshots per measurement window
nT = 100;                           % Number of windows in the first 1,000 recorded snapshots

methods = {'MUSIC','Root-MUSIC','GS-MUSIC','GS-root-MUSIC', ...
    'JWS-MUSIC','JWS-root-MUSIC'};
RMSE = measured_rmse(dataFile, algorithmDir, snapshot, nT);

figure('Position', [100,100,1000,440]);
b = bar(1:length(methods), RMSE, 0.6, 'FaceColor', 'flat', 'EdgeColor', 'none');

% Use the same color order for the bars as for the curves in main_snapshots.m.
ax = gca;
colors = ax.ColorOrder;
colorIndex = mod(0:length(methods)-1, size(colors,1)) + 1;
b.CData = colors(colorIndex,:);

xticks(1:length(methods));
xticklabels(methods);
xtickangle(0);
ylabel('RMSE (degree)');
upperLimit = ceil(2*1.15*max(RMSE))/2;
ylim([0,upperLimit]);
yticks(linspace(0,upperLimit,6));
text(1:length(methods), RMSE+0.025*max(RMSE), compose('%.3f', RMSE), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
grid on;
ax.XGrid = 'off';
end

function RMSE = measured_rmse(dataFile, algorithmDir, snapshot, nT)
DOA = [0,5];
P = numel(DOA);
Ls = 9;                             % Subarray length for forward-backward spatial smoothing (FBSS)
beta = 0.5;
roll = 1;
resample_times = 50;

oldFolder = pwd;
oldPath = path;
oldRng = rng;
cleanup = onCleanup(@() restore_environment(oldFolder, oldPath, oldRng)); %#ok<NASGU>
addpath(algorithmDir, '-begin');
cd(algorithmDir);

X = read_recorded_snapshots(dataFile, snapshot*nT);
rng(2026, 'twister');
doa = nan(P,nT,6);

for iT = 1:nT
    idx = (iT-1)*snapshot+(1:snapshot);
    Xgroup = X(:,idx);
    Y = fbss_snapshots(Xgroup, Ls);

    % Generate pseudo-random noise at the physical sensors using the smoothed noise power as a reference.
    % Apply the same FBSS transform to the noise and the measured snapshots.
    rawEigenvalues = svd(Xgroup*Xgroup'/snapshot);
    fbssEigenvalues = svd(Y*Y'/size(Y,2));
    rawNoisePower = mean(rawEigenvalues(P+1:end));
    fbssNoisePower = mean(fbssEigenvalues(P+1:end));
    assert(rawNoisePower>0, 'Window %d has zero estimated noise power.', iT);
    sensorRoll = roll*fbssNoisePower/rawNoisePower;
    sensorNoise = Generate_pseudo_noise(Xgroup, P, resample_times, sensorRoll);
    pseudo_noise = zeros(Ls,size(Y,2),resample_times,'like',Y);
    for i = 1:resample_times
        pseudo_noise(:,:,i) = fbss_snapshots(sensorNoise(:,:,i), Ls);
    end

    doa(:,iT,1) = MUSIC(Y, P);
    doa(:,iT,2) = rMUSIC(Y, P);
    doa(:,iT,3) = GS_MUSIC(Y, P, DOA, resample_times, roll, pseudo_noise);
    doa(:,iT,4) = GS_root_MUSIC(Y, P, DOA, resample_times, roll, pseudo_noise);
    doa(:,iT,5) = JWS_MUSIC(Y, P, resample_times, beta, roll, pseudo_noise);
    doa(:,iT,6) = JWS_root_MUSIC(Y, P, resample_times, beta, roll, pseudo_noise);
end

assert(all(isfinite(doa),'all'), 'At least one method returned a nonfinite DOA.');
RMSE = zeros(1,size(doa,3));
for i = 1:numel(RMSE)
    % Sum the squared DOA errors within each window and average over windows.
    RMSE(i) = sqrt(rmse(doa(:,:,i), DOA));
end
end

function X = read_recorded_snapshots(dataFile, count)
data = load(dataFile, 'x');
x = double(data.x);
x = x-mean(x,1);

fs = 48e3;
f0 = 4e3;
nfft = 480;
hop = 240;
assert(size(x,2)==16 && size(x,1)>=(count-1)*hop+nfft, ...
    'The recording must contain 16 channels and at least %d snapshots.', count);
n = (0:nfft-1).';
window = 0.5-0.5*cos(2*pi*n/(nfft-1));
bin = round(f0*nfft/fs);
tone = window.*exp(-1j*2*pi*n*bin/nfft);

% Conjugate each snapshot to match the positive phase convention of the steering vectors.
X = complex(zeros(size(x,2),count));
for i = 1:count
    idx = (i-1)*hop+(1:nfft);
    X(:,i) = conj(x(idx,:).'*tone);
end
end

function Y = fbss_snapshots(X, Ls)
[M,T] = size(X);
nSubarrays = M-Ls+1;
Y = zeros(Ls,2*nSubarrays*T,'like',X);
for i = 1:nSubarrays
    F = X(i:i+Ls-1,:);
    cols = (i-1)*2*T+(1:T);
    Y(:,cols) = F;
    Y(:,cols+T) = flipud(conj(F));
end
end

function restore_environment(oldFolder, oldPath, oldRng)
cd(oldFolder);
path(oldPath);
rng(oldRng);
end
