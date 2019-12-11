function Run(go)
% GridOp.Run
%
% Description: do the next gridop run
%
% Syntax: go.Run
%
% Updated: 2013-09-24
% Copyright 2013 Alex Schlegel (schlegel@gmail.com).  This work is licensed
% under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
% License.
nRun	= GO.Param('exp','runs');
tRun	= GO.Param('trrun');
nTrial	= GO.Param('trialperrun');

trPrePost	= GO.Param('time','prepost');
trTrial		= GO.Param('trtrial');
trFeedback	= GO.Param('time','feedback');
trRest		= GO.Param('time','rest');


%get the current run
	kRun	= go.Experiment.Info.Get('go','run');
	kRun	= ask('Next run','dialog',false,'default',kRun);

%add to the log
	go.Experiment.AddLog(['run ' num2str(kRun) ' start']); 

%add to info
	go.Experiment.Info.Set('go','run',kRun);

%perform the run
	%disable the keyboard
		ListenChar(2);
	%show a status
		go.Mapping('wait',false);
		go.Experiment.Window.Flip('waiting for scanner');
		
		go.Experiment.Scanner.StartScan(tRun);
	%let idle processes execute
		go.Experiment.Scheduler.Wait;
	%open the prompt texture
		go.Experiment.Window.OpenTexture('prompt');
	%do the run
		%set up the sequence
			cF			=	[	{@DoRest}
								repmat({@DoTrial; @DoFeedback; @DoRest},[nTrial 1])
							];
			tSequence	= 	cumsum([
								trPrePost
								repmat([trTrial; trFeedback; trRest],[nTrial-1 1])
								trTrial
								trFeedback
								trPrePost
							]) + 1;
		%execute the sequence
			kTrial		= 0;
			nCorrect	= 0;
			res			= [];
			[tStart,tEnd,tSequenceActual]	= go.Experiment.Sequence.Linear(cF,tSequence,'tstart',1,'tbase','absolute');
		%save the results
			result			= go.Experiment.Info.Get('go','result');
			result{kRun}	= [result{kRun} res];
			go.Experiment.Info.Set('go','result',result);
	%scanner stopped
		go.Experiment.Scanner.StopScan;
	%blank the screen
		go.Experiment.Show.Text('<color:red><size:3>RELAX!</size></color>');
		go.Experiment.Window.Flip;
	%close the prompt texture
		go.Experiment.Window.CloseTexture('prompt');
	%enable the keyboard
		ListenChar(1);


%add to the log
	go.Experiment.AddLog(['run ' num2str(kRun) ' end']);
%save
	go.Experiment.Info.Save;

%increment run or end
	if kRun < nRun
		go.Experiment.Info.Set('go','run',kRun+1);
	else
		if isequal(ask('End experiment?','dialog',false,'choice',{'y','n'}),'y')
			go.End;
		else
			disp('*** Remember to go.End ***');
		end
	end

%------------------------------------------------------------------------------%
function tNow = DoRest(tNow,tNext)
	go.Experiment.AddLog('rest');
	
	%blank the screen
		go.Experiment.Show.Blank;
		go.Experiment.Window.Flip;
	%prepare the next prompt
		if kTrial<nTrial
			kTrial	= kTrial + 1;
			
			go.ShowPrompt(kRun,kTrial,'window','prompt');
		end
	
	go.Experiment.Scheduler.Wait;
end
%------------------------------------------------------------------------------%
function tNow = DoTrial(tNow,tNext)
	%execute the trial
		go.Experiment.AddLog(['trial ' num2str(kTrial)]);
		
		resCur	= go.Trial(kRun,kTrial,tNow,...
					'prompttexture','prompt'	  ...
					);
		
		if isempty(res)
			res	= resCur;
		else
			res(end+1)	= resCur;
		end
end
%------------------------------------------------------------------------------%
function tNow = DoFeedback(tNow,tNext)
	%add a log message
		nCorrect	= nCorrect + res(end).correct;
		strCorrect	= conditional(res(end).correct,'y','n');
		strTally	= [num2str(nCorrect) '/' num2str(kTrial)];
	
	go.Experiment.AddLog(['feedback (' strCorrect ', ' strTally ')']);
	
	%get the message and change in winnings
		if res(end).correct
			strFeedback	= 'Yes!';
			strColor	= 'green';
			dWinning	= GO.Param('rewardpertrial');
		else
			strFeedback	= 'No!';
			strColor	= 'red';
			dWinning	= -GO.Param('penaltypertrial');
		end
	%update the winnings
		go.reward	= max(go.reward + dWinning,GO.Param('reward','base'));
	
	strText	= ['<color:' strColor '>' strFeedback ' (' StringMoney(dWinning,'sign',true) ')</color>\n\nCurrent total: ' StringMoney(go.reward)]; 
	
	go.Experiment.Show.Text(strText);
	go.Experiment.Window.Flip;
end
%------------------------------------------------------------------------------%

end
