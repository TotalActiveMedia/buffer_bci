if ( exist('initPaths','file') ) 
  initPaths;
else
  run ../utilities/initPaths;
end

buffhost='localhost';buffport=1972;
global ft_buff; ft_buff=struct('host',buffhost,'port',buffport);
% wait for the buffer to return valid header information
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

% set the real-time-clock to use
initgetwTime();
initsleepSec();
% init the buffer clock alignment
global rtclockrecord rtclockmb;
[rtclockmb rtclockrecord]=buffer_alignrtClock();
clockUpdateTime=getwTime();
clockUpdateInterval=1; %


% make the target sequence
sentences={'hello world','this is new!','BCI is fun!'};
interSentenceDuration=3;
interCharDuration=1;


% ----------------------------------------------------------------------------
%    FILL IN YOUR CODE BELOW HERE
% ----------------------------------------------------------------------------


% useful functions


% make the stimulus, i.e. put a text box in the middle of the axes
clf;
set(gcf,'color',[0 0 0],'toolbar','none','menubar','none'); % black figure
set(gca,'visible','off','color',[0 0 0]); % black axes
h=text(.5,.5,'text','HorizontalAlignment','center','VerticalAlignment','middle',...
       'FontUnits','normalized','fontsize',.2,'color',[1 1 1],'visible','off'); 
% update the text displayed
set(h,'string','new string');

% send event annotating the current time
sendEvent('stimulus.sentences','start');

% sleep (accuratly) for a certain duration
sleepSec(interCharDuration);

% keep the clock in sync  
ftime=getwTime();
if ( (ftime-clockUpdateTime)>clockUpdateInterval ) % keep the clock sync
  status=buffer('wait_dat',[-1 -1 -1]); % current sample info
  [rtclockmb rtclockrecord]=updateClocks(rtclockrecord,status.nsamples,getwTime()); 
  clockUpdateTime=ftime;
end

% wait for a key press
msg=msgbox({'Press OK to continue'},'Continue?');while ishandle(msg); pause(.2); end;
