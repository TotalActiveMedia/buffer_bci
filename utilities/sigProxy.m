function []=sigProxy(withServer,buffhost,buffport,varargin);
% simple eegviewer function
%
% eegViewer(buffhost,buffport,varargin)
%
% Inputs:
%  buffhost -- host name where the buffer is
%  buffport -- port to connect to
% Options:
%  endType -- event type which means stop viewing eeg   ('end.training')
%  trlen_ms/trlen_samp -- amount of data to plot for each channel (5000ms)
%  updateFreq -- [single] frequency to re-draw the display           (4)
%  detrend    -- [bool]  detrend the data before plotting            (1)
%  fftfilter  -- [4x1] spectral filter to apply to the data before plotting ([.1 .3 45 47])
%  downsample -- [single] frequency to downsample to before drawing display (128)
%  spatfilt   -- 'str' name of type of spatial filtering to do to the data   ('car')
%  capFile    -- [str] capFile name to get the electrode positions from      ('1010')
%  overridechnms -- [bool] flag if we use the channel names from the capFile rather than the buffer (0)
%  welch_width_ms -- [single] size in time of the welch window               (500ms) 
%                      -> defines the frequency resolution for the frequency view of the data.   
%  freqbands  -- [2x1] frequency bands to display in the freq-domain plot    (opts.fftfilter)
%  noisebands -- [2x1] frequency bands to display for the 50 Hz noise plot   ([45 47 53

page_screen_output(0);
page_output_immediately(1); % prevent buffering output

if ( nargin<1 || isempty(withServer) ) withServer=0; end;

if ( withServer )
  server = socket(AF_INET, SOCK_STREAM, 0);
  bind(server, 9001);
  server_info = struct("addr", "192.168.43.153", "port", 9001);
  connect(server, server_info);
end;

wb=which('buffer'); if ( isempty(wb) || isempty(strfind('dataAcq',wb)) ) run('../utilities/initPaths.m'); end;
opts=struct('endType','end.training','verb',1,'trlen_ms',5000,'trlen_samp',[],'updateFreq',4,'detrend',1,'fftfilter',[.1 .3 45 47],'freqbands',[],'downsample',128,'spatfilt','car','capFile','muse.txt','overridechnms',0,'welch_width_ms',500,'noisebands',[45 47 53 55],'noiseBins',[0 1],'timeOut_ms',1000);
opts=parseOpts(opts,varargin);
if ( nargin<2 || isempty(buffhost) ) buffhost='localhost'; end;
if ( nargin<3 || isempty(buffport) ) buffport=1972; end;
if ( isempty(opts.freqbands) && ~isempty(opts.fftfilter) ) opts.freqbands=opts.fftfilter; end;

% get channel info for plotting
hdr=[];
while ( isempty(hdr) || ~isstruct(hdr) || (hdr.nchans==0) ) % wait for the buffer to contain valid data
  try 
    hdr=buffer('get_hdr',[],buffhost,buffport); 
  catch
    hdr=[];
    fprintf('Invalid header info... waiting.\n');
  end;
  pause(1);
end;
capFile=opts.capFile; overridechnms=opts.overridechnms; 
if(isempty(opts.capFile)) 
  [fn,pth]=uigetfile('../utilities/*.txt','Pick cap-file'); capFile=fullfile(pth,fn);
  if ( isequal(fn,0) || isequal(pth,0) ) capFile='1010.txt'; end; % 1010 default if not selected
end
if ( ~isempty(strfind(capFile,'1010.txt')) ) overridechnms=0; else overridechnms=1; end; % force default override
di = addPosInfo(hdr.channel_names,capFile,overridechnms); % get 3d-coords
ch_pos=cat(2,di.extra.pos2d); ch_names=di.vals; % extract pos and channels names
iseeg=[di.extra.iseeg];

if ( isfield(hdr,'fSample') ) fs=hdr.fSample; else fs=hdr.fsample; end;
trlen_samp=opts.trlen_samp;
if ( isempty(trlen_samp) && ~isempty(opts.trlen_ms) ) trlen_samp=round(opts.trlen_ms*fs/1000); end;
update_samp=ceil(fs/opts.updateFreq);
trlen_samp=ceil(trlen_samp/update_samp)*update_samp;
fprintf('tr_samp = %d update_samp = %d\n',trlen_samp,update_samp);
blkIdx = trlen_samp-update_samp; % index for where new data gets inserted
if ( isempty(opts.downsample) || opts.downsample>fs ) 
    times=(-trlen_samp+1:0)./fs;
else
    times=(-ceil((trlen_samp+1)*opts.downsample/fs):0)/opts.downsample;
end
freqs=0:1000/opts.welch_width_ms:fs/2;
[ans,freqIdx(1)]=min(abs(freqs-opts.freqbands(1))); 
[ans,freqIdx(2)]=min(abs(freqs-opts.freqbands(max(end,2))));
[ans,noiseIdx(1)]=min(abs(freqs-opts.noisebands(1))); 
[ans,noiseIdx(2)]=min(abs(freqs-opts.noisebands(max(end,2))));

% make the spectral filter
filt=[]; if ( ~isempty(opts.freqbands)) filt=mkFilter(trlen_samp/2,opts.freqbands,fs/trlen_samp);end
outsz=[trlen_samp trlen_samp];if(~isempty(opts.downsample)) outsz(2)=min(outsz(2),round(trlen_samp*opts.downsample/fs)); end;
  
% recording the ERP data
rawdat    = zeros(sum(iseeg),outsz(1));
ppdat     = zeros(sum(iseeg),outsz(2));
% and the spectrogram version
[ppspect,start_samp,freqs]=spectrogram(ppdat,2,'width_ms',opts.welch_width_ms,'fs',hdr.fsample);
ppspect=ppspect(:,freqIdx(1):freqIdx(2),:); % subset to freq range of interest
start_s=-start_samp(end:-1:1)/hdr.fsample;

fprintf('Freqs : %s',sprintf('%g, ',freqs(freqIdx(1):freqIdx(2))));
                                                                               
% make popup menu for selection of TD/FD
modehdl=[]; vistype=0;

figure();
alpha_hist = beta_hist = [];
                                                                               
endTraining=false; state=[];
cursamp=hdr.nSamples;
while ( ~endTraining )  
  % wait for new data to be available
  status=buffer('wait_dat',[cursamp+update_samp inf opts.timeOut_ms],buffhost,buffport);
  if( status.nSamples < cursamp+update_samp )
    fprintf('Buffer stall detected...\n');
    pause(1);
    cursamp=status.nSamples;
  elseif ( status.nSamples > cursamp+update_samp*2 ) % missed a whole update window
    cursamp=status.nSamples - update_samp-1; % jump to the current time
  end;
  dat   =buffer('get_dat',[cursamp+1 cursamp+update_samp],buffhost,buffport);
  cursamp = cursamp+update_samp;
  
  if ( opts.verb>0 ) fprintf('.'); end;

  % shift and insert into the data buffer
  rawdat(:,1:blkIdx)=rawdat(:,update_samp+1:end);
  rawdat(:,blkIdx+1:end)=dat.buf(iseeg,:);
  
  % pre-process the data
  ppdat = rawdat;
  if ( opts.detrend ) ppdat=detrend(ppdat,2); end;
  if ( ~isempty(opts.spatfilt) ) 
    if ( strcmpi(opts.spatfilt,'car') ) ppdat=repop(ppdat,'-',mean(ppdat,1)); end
  end

  ppdat = welchpsd(ppdat,2,'width_ms',opts.welch_width_ms,'fs',hdr.fsample,'aveType','db');
  ppdat = ppdat(:,freqIdx(1):freqIdx(2));

  % update the plot
  % Get average of channels 2 until 5
  average_power = mean(ppdat(2:5,:),1);
  % Get frequency range from 8 untill 10 Hz
  mean_alpha = (mean(average_power(5:6)) * 1000) - 500;
  % Get frequency range from 14 untill 26 Hz
  mean_beta = (mean(average_power(7:13)) * 1000) - 500;
    
  alpha_hist = [alpha_hist;mean_alpha];
  alpha_hist = alpha_hist(max(1,end-100):end);
  beta_hist = [beta_hist;mean_beta];
  beta_hist = beta_hist(max(1,end-100):end);

  % red
  plot(alpha_hist, "1");
  hold on;

  % magenta
  plot(beta_hist, "4");
  hold off;
  drawnow;
  
  send_alpha = floor(mean_alpha * 4);
                                                                               
  frame = sprintf("%d\n", send_alpha);

  if ( withServer ) send(server, frame); end;
end

disconnect(server);

return;
