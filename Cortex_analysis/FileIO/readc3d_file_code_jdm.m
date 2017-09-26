
%% use this one? need VAXD? screws things up?
% 
% filepath = 'C:\Users\RatControl\Desktop\Kevin\Scripts\spineTracking\Chihuahua_8marker_yip1.c3d';
% filepath = 'C:\Users\RatControl\Desktop\Kevin\Scripts\spineTracking\Vicon2_twotap_mid12.c3d';
filepath = 'C:\Users\RatControl\Desktop\Kevin\Scripts\spineTracking\Vicon2_twotap_longrecord.c3d';

filepath = '\\OLVECZKYLABSERV\Files\Kevin\20161208\Vicon2_mid_3.c3d';
numframes = 147000;


coord = [];
fid = fopen(filepath,'r+','ieee-le');
key1 = fread(fid,1,'int8'); %key1=freadVAXD(fid,1,'int8');
key2 = fread(fid, 1, 'int8');%key2=freadVAXD(fid,1,'int8');
fseek(fid,2,'bof');
nmarkers=fread(fid,1,'int16'); %number of markers
nanalog=fread(fid,1,'int16'); %number of analog channels x #analog frames per video frame
sframe=fread(fid,1,'int16'); %# of first video frame
eframe=fread(fid,1,'int16'); %# of last video frame
nframes=eframe - sframe + 1; %number of frames
intgap=fread(fid,1,'int16'); %maximum interpolation gap allowed (in frame)
scale=fread(fid,1,'float32'); %floating-point scale factor to convert 3D-integers to ref system units
nstart=fread(fid,1,'int16'); %starting record number for 3D point and analog data
adcperframe=fread(fid,1,'int16'); %number of analog channels per video frame
fseek(fid,0,'cof');
frame_rate=fread(fid,1,'float32'); %frequency of video data
nchannel=nanalog/adcperframe; %number of analog channels
fseek(fid,512,'bof');
% 1st group data only
dat1=fread(fid,1,'int8');
dat2=fread(fid,1,'int8');
records=fread(fid,1,'int8');
proctype=fread(fid,1,'int8');
characters=fread(fid,1,'int8'); %characters in group/parametername
idnumber=fread(fid,1,'int8');
fseek(fid,(nstart-1)*512,'bof');
nframes = numframes;
nmarkers
for i=1:nframes
for j=1:nmarkers
for k=1:4
if scale < 0
coord(i,j,k)=fread(fid,1,'float32');
else
coord(i,j,k)=fread(fid,1,'int16')*scale;
end
end;
end;
for l=1:adcperframe
for m=1:nchannel
analog(i,m,l)=fread(fid,1,'int16')-2048; %zero the analog channel by subtracting 2048 bits
end;
end;
if mod(i,10000) == 0; disp(i); end;
end;
fclose(fid);
fclose 'all';

%% read frame range
filepath = 'C:\Users\RatControl\Desktop\Kevin\CortexStuff\Vicon2_twotap_longrecord.c3d';
framerange = 266484:269966;
filepath = 'C:\Users\RatControl\Desktop\Kevin\CortexStuff\Vicon2_twotap_longrecord_task.c3d';
framerange = 150712:169530;
filepath = 'C:\Users\RatControl\Desktop\Kevin\CortexStuff\Vicon2_twotap_longrecord_actxLong.c3d';
framerange = 266577:287253;

coord = [];
fid = fopen(filepath,'r+','ieee-le');
key1 = fread(fid,1,'int8'); %key1=freadVAXD(fid,1,'int8');
key2 = fread(fid, 1, 'int8');%key2=freadVAXD(fid,1,'int8');
fseek(fid,2,'bof');
nmarkers=fread(fid,1,'int16'); %number of markers
nanalog=fread(fid,1,'int16'); %number of analog channels x #analog frames per video frame
sframe=fread(fid,1,'int16'); %# of first video frame
eframe=fread(fid,1,'int16'); %# of last video frame
nframes=eframe - sframe + 1; %number of frames
intgap=fread(fid,1,'int16'); %maximum interpolation gap allowed (in frame)
scale=fread(fid,1,'float32'); %floating-point scale factor to convert 3D-integers to ref system units
nstart=fread(fid,1,'int16'); %starting record number for 3D point and analog data
adcperframe=fread(fid,1,'int16'); %number of analog channels per video frame
fseek(fid,0,'cof');
frame_rate=fread(fid,1,'float32'); %frequency of video data
nchannel=nanalog/adcperframe; %number of analog channels
fseek(fid,512,'bof');
% 1st group data only
dat1=fread(fid,1,'int8');
dat2=fread(fid,1,'int8');
records=fread(fid,1,'int8');
proctype=fread(fid,1,'int8');
characters=fread(fid,1,'int8'); %characters in group/parametername
idnumber=fread(fid,1,'int8');
fseek(fid,(nstart-1)*512,'bof');

bytesperframe = nmarkers * 4;
fseek(fid,(nstart-1)*512 + bytesperframe*framerange(1),'bof');

%** above fails, temp try this
for i = 1:framerange(1)
    for j = 1:nmarkers
        for k = 1:4
            if scale < 0; 
                fread(fid,1,'float32');
            else
                fread(fid,1,'int16')*scale;
            end
        end
    end
    if mod(i,1000)==0; disp(i); end;
end
%**

for i=1:length(framerange)%1:nframes
for j=1:nmarkers
for k=1:4
if scale < 0
coord(i,j,k)=fread(fid,1,'float32');
else
coord(i,j,k)=fread(fid,1,'int16')*scale;
end
end;
end;
for l=1:adcperframe
for m=1:nchannel
analog(i,m,l)=fread(fid,1,'int16')-2048; %zero the analog channel by subtracting 2048 bits
end;
end;
end;
fclose(fid);
fclose 'all';

%% plot 3d data

%figure; hold on;
%scatter3(coord(1,:,1),coord(1,:,2),coord(1,:,3));
%coord1 = coord(1:10000,:,1:3);
coord1 = coord;
figure; hold on;
set(gca,'Xlim',[-150 300],'Ylim', [-300 300], 'Zlim', [0 300]);
for j = 1:length(coord1);
    scatter3(coord1(j,:,1),coord1(j,:,2),coord1(j,:,3));
    %set(gca,'Xlim',[-75 125],'Ylim', [0 220], 'Zlim', [0 100]);
    set(gca,'Xlim',[-150 300],'Ylim', [-300 300], 'Zlim', [0 300]);
    title(['frame ' num2str(j)]);
    pause(.001);
    clf;
end

%% run in bparhmm
