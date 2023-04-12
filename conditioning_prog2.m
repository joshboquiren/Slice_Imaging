%% Runs session
% This function is called upon clicking the Start button in the GUI
% Requires variables param and s (serial object).
running = true;     
n=timer 

cla;    % Clear the axes before starting
%% Parameters

truncITI = min(maxITI,3*minITI);                     % minITI is really the mean ITI
licksinit = ceil(sessionLim*(truncITI+t_fxd)*10/1E3);  % number of licks to initialize = number of trials*max time per trial in s*10Hz (of licking)
cuesinit = sessionLim;                               % number of cues to initialize
cuesminusinit = numCSminus;                          % number of CS- cues to initialize
logInit = 10^6;                                      % Log of all serial input events from the Arduino
bgdpumpsinit = ceil(sessionLim*truncITI*3/T_bgd);      % number of background pumps to initialize = total time spent in ITI*rate of rewards*3. It won't be more than 3 times the expected rate

xWindow = [-(truncITI+1000) t_fxd+1000];  % Defines x-axis limits for the plot. Stretches between truncITI+1 s before cue to t_fxd+1 s after cue
fractionPSTHdisplay = 0.15;             % What fraction of the display is the PSTH?
yOffset = ceil(fractionPSTHdisplay*sessionLim/(1-fractionPSTHdisplay));% amount by which the y-axis extends beyond trials so as to plot the PSTH of licks
binSize = 1000;                         % in ms
xbins = xWindow(1):binSize:xWindow(2);  % Bins in x-axis for PSTH

ticks = -(truncITI+1000):10000:t_fxd+1000;% tick marks for x-axis of raster plot. moves through by 2s
labels = ticks'/1000;                     % convert tick labels to seconds
labelsStr = cellstr(num2str(labels));     % convert to cell of strings

durationtrialpartitionnocues = 20E3;      % When nocuesflag=1, how long should a single row for raster plot be?
%% Prep work

% Inputs from Arduino along with their "code" (defined below)
% 1 = Lick onset
% 2 = Lick offset
% 3 = Background pump
% 4 = Fixed pump
% 5 = CS+ Cue
% 6 = CS- Cue
% 7 = laser onset
% 21 = left lever press
% 22 = right lever press
% 23= trial time limit reached
% 50 = manual reward delivered
% 80 = start of new trial



% initialize arrays for licks and cues
licks = zeros(licksinit,1);
r = 0;% Counter for licks
lickoffsets = zeros(licksinit,1);
roff = 0;% Counter for lick offsets
bgdpumps = zeros(bgdpumpsinit,1);
bgdus = 0;% Counter for background rewards
fxdpumps = zeros(cuesinit,1);
fxdus = 0;% Counter for fixed rewards
cues = zeros(cuesinit,1);
cs = 0;% Counter for cues
cuesminus = zeros(cuesminusinit,1);
csminus = 0;% Counter for CS- cues
eventlog = zeros(logInit,4);
l = 0;% Counter for logged events

% initialize arrays for operant conditioning responses and DS. Added by
% ITP, dec 2016
trialPressLeft= zeros(trialNumLim,1);
trialPressRight= zeros(trialNumLim,1);
TrialNum=1;
manRWD = 0;

% The following variables are declared because they are stored in the
% buffer until a trial ends. This is done so that the plot can be aligned
% to the cue for all trials. Real time plotting cannot work in this case
% since there is variability in intercue interval.
templicks = NaN(ceil(licksinit/sessionLim),1);
templicksct = 0;                                                % count of temp licks. Calculated explicitly to speed up indexing
templicksPSTHplus = NaN(ceil(licksinit/(sessionLim-numCSminus)),sessionLim-numCSminus); % Array in which all CS+ licks are stored for calculating PSTH 
if numCSminus~=0
    templicksPSTHminus = NaN(ceil(licksinit/numCSminus),numCSminus); % Array in which all CS- licks are stored for calculating PSTH 
elseif numCSminus==0
    templicksPSTHminus = [];
end
hPSTHplus = [];                                                 % Handle to PSTH plot on CS+ trials
hPSTHminus = [];                                                % Handle to PSTH plot on CS- trials
temprewards = NaN(ceil((bgdpumpsinit+cuesinit)/sessionLim),1);
temprewardsct = 0;
tempcueplus = NaN(1,1);
tempcueminus = NaN(1,1);

if nocuesflag == 1
    templicks = [];
    temprewards = [];
end

%%% run if operant experiment, added by ITP dec 2016
if OperantExperiment==1
    
    %% Load to arduino

startT = clock;                                     % find time of start
startTStr = sprintf('%d:%d:%02.0f', ...
                    startT(4),startT(5),startT(6)); % format time
set(handles.startText,'String',startTStr)           % display time
drawnow

wID = 'MATLAB:serial:fscanf:unsuccessfulRead';      % warning id for serial read timeout
warning('off',wID)                                  % suppress warning

                                % variable to control program
fprintf('Executing trial %d\n',1);
%

try

%% Collect data from arduino
while running

  
    %set(handles.CTrialTime, 'String', num2str(time))
 
    read = [];
        if s.BytesAvailable > 0 % is data available to read? This avoids the timeout problem
            read = fscanf(s,'%u');% scan for data sent only when data is available
        display(read);
        end
        if isempty(read)
            drawnow
            continue
        end
        code = read(1);                             % read identifier for data
                
        if code == 0                                % signifies "end of session"
            break
        end
        
        l = l + 1;
        eventlog(l,:) = read;
        
        time = read(2);                             % record timestamp
        
        
      
        
        norewardflag = read(3);                     % if =1, no reward was actually given. Indicates CS- trial
        % record data as either CS+ cue, CS- cue, lick, fixed reward or
        % background reward
        if code == 1                                % Lick onset; BLACK
            if nocuesflag == 0                      % Store lick timestamp for later plotting after trial ends
                r = r + 1;
                set(handles.licksEdit,'String',num2str(r))
                licks(r) = time;                        % record timestamp of lick
                templicksct = templicksct+1;
                templicks(templicksct) = time;
            elseif nocuesflag == 1 %If only Poisson rewards are given, plot when lick occurs in real time
                r = r + 1;
                set(handles.licksEdit,'String',num2str(r))
                licks(r) = time;                        % record timestamp of lick
                trial = floor(time/durationtrialpartitionnocues);
                temptrialdur = trial*durationtrialpartitionnocues;                
             end
        elseif code == 2                            % Lick offset
            roff = roff + 1;
            lickoffsets(r) = time;
        elseif code == 3                            % Background reward; BLUE
            if nocuesflag == 0
                bgdus = bgdus + 1;            
                bgdpumps(bgdus,1) = time;
                set(handles.bgdpumpsEdit,'String',num2str(bgdus))
                temprewardsct = temprewardsct+1;
                temprewards(temprewardsct) = time;
            elseif nocuesflag == 1 %If only Poisson rewards are given, plot when reward occurs
                bgdus = bgdus + 1;            
                bgdpumps(bgdus,1) = time;
                set(handles.bgdpumpsEdit,'String',num2str(bgdus))
                trial = floor(time/durationtrialpartitionnocues);
                temptrialdur = trial*durationtrialpartitionnocues;                
             end
        elseif code == 4                            % Fixed reward; BLUE
            if norewardflag == 0                   % Indicates reward was delivered
                TrialNum = TrialNum+1;
                fprintf('Executing trial %d\n',TrialNum);                        
            set(handles.TrialNum, 'String',num2str(TrialNum))
            trialPressLeft = 0;
            trialPressRight = 0;
            set(handles.trialPressLeft,'String', num2str(trialPressLeft))
            set(handles.trailPressRight,'String', num2str(trialPressRight))
            drawnow
                fxdus = fxdus + 1;            
                fxdpumps(fxdus,1) = time;
                set(handles.fxdpumpsEdit,'String',num2str(fxdus))
                temprewardsct = temprewardsct+1;
                temprewards(temprewardsct) = time;
            end
            
            
            
            
        elseif code == 21 %if left lever was pressed
            if norewardflag==1 %if no reward was given
                trialPressLeft= trialPressLeft+1;
            set(handles.trialPressLeft,'String', num2str(trialPressLeft))
            drawnow
            end
            
        elseif code == 22 %if right lever was pressed
            if norewardflag==1 %if no reward was given
                trialPressRight= trialPressRight+1;
            set(handles.trailPressRight,'String', num2str(trialPressRight))
            drawnow
            end
            
            
        elseif code == 23
            TrialNum=TrialNum+1;
            % fprintf('Executing trial %d\n',TrialNum);
            set(handles.trialPressLeft,'String','0')
            set(handles.trailPressRight, 'String', '0')
            set(handles.TrialNum, 'String',num2str(TrialNum))
            trialPressLeft=0;
            trialPressRight=0;
            drawnow
            
            elseif code == 80            
            set(handles.TrialNum, 'String',num2str(TrialNum))
            set(handles.trialPressLeft,'String','0')
            set(handles.trailPressRight, 'String', '0')
                        trialPressLeft=0;
            trialPressRight=0;
            
            
        elseif code == 50
            manRWD= manRWD+1;
            set(handles.bgdpumpsEdit, 'String', num2str(manRWD))
            % TrialNum=TrialNum+1;
            %fprintf('Executing trial %d\n',TrialNum);
            set(handles.trialPressLeft,'String','0')
            set(handles.trailPressRight, 'String', '0')
            set(handles.TrialNum, 'String',num2str(TrialNum))
            drawnow
        
            
            temprewards = temprewards-time+t_fxd; %find timestamps wrt the cue
            tempcueplus = tempcueplus-time+t_fxd;
            tempcueminus = tempcueminus-time+t_fxd;
            templicks = templicks-time+t_fxd;
            if ~isnan(tempcueplus)
                templicksPSTHplus(1:length(templicks),cs-csminus) = templicks;
            elseif ~isnan(tempcueminus)
                templicksPSTHminus(1:length(templicks),csminus) = templicks;
            end
            
            % Raster plot
            plot([templicks templicks],[-(cs-1) -cs],'k','LineWidth',1);hold on
            plot([temprewards temprewards],[-(cs-1) -cs],'c','LineWidth',2);hold on
            plot([tempcueplus tempcueplus],[-(cs-1) -cs],'g','LineWidth',2);hold on
            plot([tempcueminus tempcueminus],[-(cs-1) -cs],'r','LineWidth',2);hold on
            
            % PSTH
            delete(hPSTHplus);delete(hPSTHminus); %Clear previous PSTH plots            
            nPSTHplus = histc(templicksPSTHplus(~isnan(templicksPSTHplus)),xbins); % Count licks in each bin for all trials until now
            nPSTHplus = nPSTHplus/max(nPSTHplus); % Plot PSTH for CS+ scaled to the available range on the y-axis            
%             hPSTHplus = bar(xbins,nPSTHplus*yOffset,1,'FaceColor',[0.47 0.67 0.19],'EdgeColor','k'); %plot histogram of lick density on CS+
            hPSTHplus = plot(xbins,nPSTHplus*yOffset,'Marker','o','MarkerFaceColor',[0.47 0.67 0.19],'Color',[0.47 0.67 0.19]);
            hold on;
            if ~isempty(templicksPSTHminus)
                nPSTHminus = histc(templicksPSTHminus(~isnan(templicksPSTHminus)),xbins); % Count licks in each bin for all trials until now
                nPSTHminus = nPSTHminus/max(nPSTHminus); % Plot PSTH for CS- scaled to the available range on the y-axis
%                 hPSTHminus = bar(xbins,nPSTHminus*yOffset,1,'FaceColor',[1 0.6 0.78],'EdgeColor','k'); %plot histogram of lick density on CS-
                hPSTHminus = plot(xbins,nPSTHminus*yOffset,'Marker','o','MarkerFaceColor',[1 0.6 0.78],'Color',[1 0.6 0.78]);
            end
            drawnow
            
            % Re-initialize the temp variables
            templicks = NaN(ceil(licksinit/sessionLim),1);
            templicksct = 0; %count of temp licks. Calculated explicitly to speed up indexing
            temprewards = NaN(ceil((bgdpumpsinit+cuesinit)/sessionLim),1);
            temprewardsct = 0;
            tempcueplus = NaN(1,1);
            tempcueminus = NaN(1,1);
        elseif code == 5                            % CS+ cue onset; GREEN
            if  norewardflag == 0 %if reward was given
                %TrialNum = TrialNum+1;
                %fprintf('Executing trial %d\n',TrialNum);                        
            set(handles.TrialNum, 'String',num2str(TrialNum))
            trialPressLeft = 0;
            set(handles.trialPressLeft,'String', num2str(trialPressLeft))
            drawnow
            end
            cs = cs + 1;            
            cues(cs+1,1) = time;                      % record timestamp of tone
            set(handles.cuesEdit,'String',num2str(cs-csminus))
            tempcueplus = time;
            if cs<sessionLim
                %fprintf('Executing trial %d\n',TrialNum);
            end
        elseif code == 6                            % CS- cue onset; RED
            cs = cs + 1;            
            cues(cs+1,1) = time;                      % record timestamp of tone
            TrialNum = TrialNum+1;
            
            csminus = csminus + 1;            
            cuesminus(csminus+1,1) = time;            % record timestamp of tone
            set(handles.cuesminusEdit,'String',num2str(csminus))            
            tempcueminus = time;
            fprintf('Executing trial %d\n',TrialNum);
            
            if cs<sessionLim
                %fprintf('Executing trial %d\n',cs+1);
            end    
   
        end 
       
end 


       

 %% Save operant data

    format = 'yymmdd-HHMMSS';
    date = datestr(now,format);
    
    %if nocuesflag == 1
        str = 'Poisson_';
    %elseif nocuesflag ==0
        str = [];
    %end
    if csminusprob==100 && csplusprob==0
        revstr = 'reverse_';
    else
        revstr = [];
    end
    
    %if lasercueflag == 1 && laserrewardflag == 0
        laserstr = 'lasercue_';
    %elseif lasercueflag == 0 && laserrewardflag == 1
        laserstr = 'laserreward_';
    %elseif lasercueflag == 1 && laserrewardflag == 1
        laserstr = 'lasercuereward_';
    %elseif randlaserflag ==1
        laserstr = 'randlaser_';
    %else
        laserstr = [];
    %end
    
    %if trialbytrialbgdpumpflag == 1
        bgdpumpstr = 'trialbytrialbgd_';
    %else
        bgdpumpstr = [];
    %end
    
    if r_fxd == 0
        extinctionstr = 'extinction_';
    else
        extinctionstr = [];
    end
    
    if csplusprob+csminusprob ~= 100
        probstr = ['plusprob' num2str(csplusprob) 'minusprob' num2str(csminusprob) '_'];
    else
        probstr = [];
    end
    
    assignin('base','eventlog',eventlog);
    file = [saveDir fname '_' num2str(r_bgd) '_' num2str(T_bgd) '_'  str revstr probstr laserstr bgdpumpstr extinctionstr date '.mat'];
    
    save(file, 'licks', 'lickoffsets', 'cues', 'cuesminus', 'bgdpumps', 'fxdpumps', 'eventlog', 'params')
catch exception
    if r < licksinit
        licks = licks(1:r);
    end
    if roff < licksinit
        lickoffsets = lickoffsets(1:roff);
    end
    if fxdus < cuesinit
        fxdpumps = fxdpumps(1:fxdus,:);
    end
    if bgdus < bgdpumpsinit
        bgdpumps = bgdpumps(1:bgdus,:);
    end
    if cs < cuesinit
        cues = cues(1:cs,:);
    end
    if csminus < cuesminusinit
        cuesminus = cuesminus(1:csminus,1);
    end
    if l < logInit
        eventlog = eventlog(1:l,:);
    end    
end
end


%% run if classical conditioning experiment (VJ's original code), conditioonal added by ITP dec 2016

if OperantExperiment==0

% setup plot
axes(actvAx)                            % make the activity axes the current one
if nocuesflag == 0
    plot(xWindow,[0 0],'k','LineWidth',2);hold on                   % start figure for plots
    set(actvAx,'ytick',[], ...
               'ylim',[-sessionLim yOffset+1], ...
               'ytick',[], ...
               'xlim',xWindow, ...
               'xtick',ticks, ...
               'xticklabel',labelsStr');        % set labels: Raster plot with y-axis containing trials. Chronological order = going from top to bottom
    xlabel('time (s)');
    ylabel('Trials');
elseif nocuesflag == 1
    plot([0 0;0 0],[0 0;-1 -1],'w');hold on
    xlabel('time (s)');
    ylabel(' ');
    xlim([-1000 durationtrialpartitionnocues+1000]);
    set(actvAx,'ytick',[],...
               'xtick',0:2000:durationtrialpartitionnocues,...
               'XTickLabel',num2str((0:2000:durationtrialpartitionnocues)'/1000));
end

drawnow


%% Load to arduino

startT = clock;                                     % find time of start
startTStr = sprintf('%d:%d:%02.0f', ...
                    startT(4),startT(5),startT(6)); % format time
set(handles.startText,'String',startTStr)           % display time
drawnow

wID = 'MATLAB:serial:fscanf:unsuccessfulRead';      % warning id for serial read timeout
warning('off',wID)                                  % suppress warning

running = true;                                     % variable to control program
fprintf('Executing trial %d\n',1);
while running 
    if (n1>1)
        continue 
    end 
%%
try
    
%% Collect data from arduino
    while running
        read = [];
        if s.BytesAvailable > 0 % is data available to read? This avoids the timeout problem
            read = fscanf(s,'%u');% scan for data sent only when data is available
            
        end
        if isempty(read)
            drawnow
            continue
        end
        
        
        l = l + 1;
        eventlog(l,:) = read;
        
        time = read(2);                             % record timestamp
        
        norewardflag = read(3);                     % if =1, no reward was actually given. Indicates CS- trial

        % record data as either CS+ cue, CS- cue, lick, fixed reward or
        % background reward
        if code == 1                                % Lick onset; BLACK
            if nocuesflag == 0                      % Store lick timestamp for later plotting after trial ends
                r = r + 1;
                set(handles.licksEdit,'String',num2str(r))
                licks(r) = time;                        % record timestamp of lick
                templicksct = templicksct+1;
                templicks(templicksct) = time;
            elseif nocuesflag == 1 %If only Poisson rewards are given, plot when lick occurs in real time
                r = r + 1;
                set(handles.licksEdit,'String',num2str(r))
                licks(r) = time;                        % record timestamp of lick
                trial = floor(time/durationtrialpartitionnocues);
                temptrialdur = trial*durationtrialpartitionnocues;                
                plot([time-temptrialdur;time-temptrialdur],[-trial;-trial-1],'k','LineWidth',1);hold on
            end
        elseif code == 2                            % Lick offset
            roff = roff + 1;
            lickoffsets(r) = time;
        elseif code == 3                            % Background reward; BLUE
            if nocuesflag == 0
                bgdus = bgdus + 1;            
                bgdpumps(bgdus,1) = time;
                set(handles.bgdpumpsEdit,'String',num2str(bgdus))
                temprewardsct = temprewardsct+1;
                temprewards(temprewardsct) = time;
            elseif nocuesflag == 1 %If only Poisson rewards are given, plot when reward occurs
                bgdus = bgdus + 1;            
                bgdpumps(bgdus,1) = time;
                set(handles.bgdpumpsEdit,'String',num2str(bgdus))
                trial = floor(time/durationtrialpartitionnocues);
                temptrialdur = trial*durationtrialpartitionnocues;                
                plot([time-temptrialdur;time-temptrialdur],[-trial;-trial-1],'c','LineWidth',2);hold on
            end
        elseif code == 4                            % Fixed reward; BLUE
            if norewardflag == 0                      % Indicates CS+ trial
                fxdus = fxdus + 1;            
                fxdpumps(fxdus,1) = time;
                set(handles.fxdpumpsEdit,'String',num2str(fxdus))
                temprewardsct = temprewardsct+1;
                temprewards(temprewardsct) = time;
            end
            
            temprewards = temprewards-time+t_fxd; %find timestamps wrt the cue
            tempcueplus = tempcueplus-time+t_fxd;
            tempcueminus = tempcueminus-time+t_fxd;
            templicks = templicks-time+t_fxd;
            if ~isnan(tempcueplus)
                templicksPSTHplus(1:length(templicks),cs-csminus) = templicks;
            elseif ~isnan(tempcueminus)
                templicksPSTHminus(1:length(templicks),csminus) = templicks;
            end
            
            % Raster plot
            plot([templicks templicks],[-(cs-1) -cs],'k','LineWidth',1);hold on
            plot([temprewards temprewards],[-(cs-1) -cs],'c','LineWidth',2);hold on
            plot([tempcueplus tempcueplus],[-(cs-1) -cs],'g','LineWidth',2);hold on
            plot([tempcueminus tempcueminus],[-(cs-1) -cs],'r','LineWidth',2);hold on
            
            % PSTH
            delete(hPSTHplus);delete(hPSTHminus); %Clear previous PSTH plots            
            nPSTHplus = histc(templicksPSTHplus(~isnan(templicksPSTHplus)),xbins); % Count licks in each bin for all trials until now
            nPSTHplus = nPSTHplus/max(nPSTHplus); % Plot PSTH for CS+ scaled to the available range on the y-axis            
%             hPSTHplus = bar(xbins,nPSTHplus*yOffset,1,'FaceColor',[0.47 0.67 0.19],'EdgeColor','k'); %plot histogram of lick density on CS+
            hPSTHplus = plot(xbins,nPSTHplus*yOffset,'Marker','o','MarkerFaceColor',[0.47 0.67 0.19],'Color',[0.47 0.67 0.19]);
            hold on;
            if ~isempty(templicksPSTHminus)
                nPSTHminus = histc(templicksPSTHminus(~isnan(templicksPSTHminus)),xbins); % Count licks in each bin for all trials until now
                nPSTHminus = nPSTHminus/max(nPSTHminus); % Plot PSTH for CS- scaled to the available range on the y-axis
%                 hPSTHminus = bar(xbins,nPSTHminus*yOffset,1,'FaceColor',[1 0.6 0.78],'EdgeColor','k'); %plot histogram of lick density on CS-
                hPSTHminus = plot(xbins,nPSTHminus*yOffset,'Marker','o','MarkerFaceColor',[1 0.6 0.78],'Color',[1 0.6 0.78]);
            end
            drawnow

            % Re-initialize the temp variables
            templicks = NaN(ceil(licksinit/sessionLim),1);
            templicksct = 0; %count of temp licks. Calculated explicitly to speed up indexing
            temprewards = NaN(ceil((bgdpumpsinit+cuesinit)/sessionLim),1);
            temprewardsct = 0;
            tempcueplus = NaN(1,1);
            tempcueminus = NaN(1,1);
        elseif code == 5                            % CS+ cue onset; GREEN
            cs = cs + 1;            
            cues(cs+1,1) = time;                      % record timestamp of tone
            set(handles.cuesEdit,'String',num2str(cs-csminus))
            tempcueplus = time;
            if cs<sessionLim
                %fprintf('Executing trial %d\n',cs+1);
            end
        elseif code == 6                            % CS- cue onset; RED
            cs = cs + 1;            
            cues(cs+1,1) = time;                      % record timestamp of tone
            
            csminus = csminus + 1;            
            cuesminus(csminus+1,1) = time;            % record timestamp of tone
            set(handles.cuesminusEdit,'String',num2str(csminus))            
            tempcueminus = time;
            
            if cs<sessionLim
                %fprintf('Executing trial %d\n',cs+1);
            end    
        end
    end

    if r < licksinit
        licks = licks(1:r);
    end
    if roff < licksinit
        lickoffsets = lickoffsets(1:roff);
    end
    if fxdus < cuesinit
        fxdpumps = fxdpumps(1:fxdus,:);
    end
    if bgdus < bgdpumpsinit
        bgdpumps = bgdpumps(1:bgdus,:);
    end
    if cs < cuesinit
        cues = cues(1:cs,:);
    end
    if csminus < cuesminusinit
        cuesminus = cuesminus(1:csminus,1);
    end
    if l < logInit
        eventlog = eventlog(1:l,:);
    end

%% Save data

    format = 'yymmdd-HHMMSS';
    date = datestr(now,format);
    
    if nocuesflag == 1
        str = 'Poisson_';
    elseif nocuesflag ==0
        str = [];
    end
    if csminusprob==100 && csplusprob==0
        revstr = 'reverse_';
    else
        revstr = [];
    end
    
    if lasercueflag == 1 && laserrewardflag == 0
        laserstr = 'lasercue_';
    elseif lasercueflag == 0 && laserrewardflag == 1
        laserstr = 'laserreward_';
    elseif lasercueflag == 1 && laserrewardflag == 1
        laserstr = 'lasercuereward_';
    elseif randlaserflag ==1
        laserstr = 'randlaser_';
    else
        laserstr = [];
    end
    
    if trialbytrialbgdpumpflag == 1
        bgdpumpstr = 'trialbytrialbgd_';
    else
        bgdpumpstr = [];
    end
    
    if r_fxd == 0
        extinctionstr = 'extinction_';
    else
        extinctionstr = [];
    end
    
    if csplusprob+csminusprob ~= 100
        probstr = ['plusprob' num2str(csplusprob) 'minusprob' num2str(csminusprob) '_'];
    else
        probstr = [];
    end
    
    assignin('base','eventlog',eventlog);
    file = [saveDir fname '_' num2str(r_bgd) '_' num2str(T_bgd) '_'  str revstr probstr laserstr bgdpumpstr extinctionstr date '.mat'];
    save(file, 'licks', 'lickoffsets', 'cues', 'cuesminus', 'bgdpumps', 'fxdpumps', 'eventlog', 'params')

catch exception
    if r < licksinit
        licks = licks(1:r);
    end
    if roff < licksinit
        lickoffsets = lickoffsets(1:roff);
    end
    if cs < cuesinit
        cues = cues(1:cs,:);
    end
    if csminus < cuesminusinit
        cuesminus = cuesminus(1:csminus,1);
    end
    if fxdus < cuesinit
        fxdpumps = fxdpumps(1:fxdus,:);
    end
    if bgdus < bgdpumpsinit
        bgdpumps = bgdpumps(1:bgdus,:);
    end
    if l < logInit
        eventlog = eventlog(1:l,:);
    end
    
    fprintf(s,'1');                                  % send stop signal to arduino; 49 in Arduino is the ASCII code for 1
    display('Error running program.')
    format = 'yymmdd-HHMMSS';
    date = datestr(now,format);
    
    
    if nocuesflag == 1
        str = 'Poisson_';
    elseif nocuesflag ==0
        str = [];
    end
    
    if csminusprob==100 && csplusprob==0
        revstr = 'reverse_';
    else
        revstr = [];
    end
    
    if lasercueflag == 1 && laserrewardflag == 0
        laserstr = 'lasercue_';
    elseif lasercueflag == 0 && laserrewardflag == 1
        laserstr = 'laserreward_';
    elseif lasercueflag == 1 && laserrewardflag == 1
        laserstr = 'lasercuereward_';
    elseif randlaserflag ==1
        laserstr = 'randlaser_';
    else
        laserstr = [];
    end
    
    if trialbytrialbgdpumpflag == 1
        bgdpumpstr = 'trialbytrialbgd_';
    else
        bgdpumpstr = [];
    end
    
    if r_fxd == 0
        extinctionstr = 'extinction_';
    else
        extinctionstr = [];
    end
    
    if csplusprob+csminusprob ~= 100
        probstr = ['plusprob' num2str(csplusprob) 'minusprob' num2str(csminusprob) '_'];
    else
        probstr = [];
    end
    
    assignin('base','eventlog',eventlog);
    
    file = [saveDir fname '_' num2str(r_bgd) '_' num2str(T_bgd) '_'  str revstr probstr laserstr bgdpumpstr extinctionstr date '.mat'];

    save(file, 'licks', 'lickoffsets', 'cues', 'cuesminus', 'bgdpumps', 'fxdpumps', 'eventlog', 'params','exception')
end
end
end

