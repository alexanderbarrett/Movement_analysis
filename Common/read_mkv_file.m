function read_mkv_file(cameradirectory,camerabase)
% LCamdir = 'CameraL';
% RCamdir = 'CameraR';
% LCamdir = 'Cam1';
% RCamdir = 'Cam2';
ext = '.mkv';
convext = '.mp4';

tempfolder = camerabase;

% datafolder = 'X:\Data\';
% savefolder = 'H:\Analysis\';
% datafolder = 'X:\Data\';
% savefolder = 'H:\Analysis\';
% ratName = behStruct(1).ratName;
% datadir = [datafolder ratName '\'];

ffmpegPath = '"C:\Users\Jesse Marshall\Documents\ffmpeg\bin\ffmpeg.exe"';
[ffmpegPath_dir,~,~] = fileparts(ffmpegPath);
ffmpegPath_dir = strrep(ffmpegPath_dir,'"','');


%
% % Get cam directory folder list
% Lfolddir = dir([cameradirectory '\']);
% Lfoldlist = int64(zeros(1, numel(Lfolddir)-2));
% for f = 1 : length(Lfoldlist)
%     Lfoldlist(f) = int64(str2num(['int64(' Lfolddir(f+2).name ')']));
% end
% %
% Rfolddir = dir([datadir RCamdir '\']);
% Rfoldlist = int64(zeros(1, numel(Rfolddir)-2));
% for f = 1 : length(Rfoldlist)
%     Rfoldlist(f) = int64(str2num(['int64(' Rfolddir(f+2).name ')']));
% end

% Delete temp files if they exist
% if exist([tempfolder 'Lvid' convext], 'file') == 2
%     delete([tempfolder 'Lvid' convext]);
% end

% if exist([tempfolder 'Rvid' convext], 'file') == 2
%     delete([tempfolder 'Rvid' convext]);
% end

Lflist = dir([cameradirectory '\*' ext]);
LFileFr1 = zeros(length(Lflist), 1);
for f = 1 : length(Lflist)
    LFileFr1(f) = str2double(Lflist(f).name(1:end-length(ext)));
end
LFileFr1 = sort(LFileFr1);

for ll = 1:numel(LFileFr1)
    outputfile = [tempfolder '\' num2str(LFileFr1(ll))  convext];
    if (~exist(outputfile,'file'))
        [~, ~] = system([ffmpegPath ' -i ' cameradirectory '\' num2str(LFileFr1(ll)) ext ' -vcodec copy ' tempfolder '\'  num2str(LFileFr1(ll)) convext]);
    end
end
%Lexist = exist(strcat(ffmpegPath_dir,'\Temp\',num2str(LFileFr1(1)),'.mp4'),'file');
%datastore=

%  [~, ~] = system([ffmpegPath ' -i ' datadir RCamdir '\' Rfold '\' num2str(RFileFr1(1)) ext ' -vcodec copy ' tempfolder 'Rvid' convext]);

%             LVideo = zeros(
%
%              for chunk = 1:ceil(LEDframesPerFile/readChunk)
%                     vstart = readChunk * (chunk-1) + 1;
%                     vend = min([LEDframesPerFile vstart+readChunk-1]);
%                     Lframes = read(Lvp, [vstart vend]);
%                     Lframes = double(imrotate(Lframes,-90))/255;
%                     Lvideo(:,:,:,vstart:vend) = Lframes(LledROI(2):LledROI(2)+LledROI(4),LledROI(1):LledROI(1)+LledROI(3),:,:);
%                 end
%   Rexist = exist(strcat(ffmpegPath_dir,'\Temp\','Rvid','.mp4'),'file');

end
