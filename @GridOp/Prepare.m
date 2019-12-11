function Prepare(go,varargin)
% GridOp.Prepare
%
% Description: prepare to run a gridop experiment
%
% Syntax: go.Prepare()
%
% Updated: 2013-09-18
% Copyright 2013 Alex Schlegel (schlegel@gmail.com).  This work is licensed
% under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
% License.
global strDirBase;

%stimulus and operation mapping
	bMappingExists	= ~isempty(go.Experiment.Subject.Get('map_stim'));
	if bMappingExists
		mapStim	= go.Experiment.Subject.Get('map_stim');
		mapOp	= go.Experiment.Subject.Get('map_op');
	else
		bExist	= false;
		if go.Experiment.Info.Get('experiment','debug')<2
			mapStim	= randomize((1:4)');
			mapOp	= randomize((1:4)');
			
			go.Experiment.Subject.Set('map_stim',mapStim);
			go.Experiment.Subject.Set('map_op',mapOp);
		else
			[mapStim,mapOp]	= deal({'default'});
		end 
	end
	
	strExist	= [' (' conditional(bMappingExists,'existing','new') ')'];
	go.Experiment.AddLog(['stimulus mapping: ' join(mapStim,' ') strExist]);
	go.Experiment.AddLog(['operation mapping: ' join(mapOp,' ') strExist]);
%response buttons
	kCorrect	= GO.Param('response','correct');
	kIncorrect	= GO.Param('response','incorrect');
	go.Experiment.Input.Set('response',{kCorrect,kIncorrect});
	go.Experiment.Input.Set('correct',kCorrect);
	go.Experiment.Input.Set('incorrect',kIncorrect);
%run
	go.Experiment.Info.Set('go','run',1);
	go.Experiment.Info.Set('go','result',cell(GO.Param('exp','runs'),1));
%trial info
	nRun	= GO.Param('exp','runs');
	nRep	= GO.Param('exp','reps');
	
	%condition order
		[kS,kO] 	= ndgrid(1:4,1:4);
		kCondition	= kS + 10*kO;
		
		nCondition	= numel(kCondition);
		nTrialPer	= nCondition*nRep;
		
		order	= GenOrder(kCondition);
		
		go.Experiment.Info.Set('go',{'trial','input'},decget(order,0));
		go.Experiment.Info.Set('go',{'trial','op'},decget(order,1));
	%prompt location
		go.Experiment.Info.Set('go',{'trial','prompt','location'},GenOrder(1:4));
	%target test location
		go.Experiment.Info.Set('go',{'trial','test','location'},GenOrder(1:4));
	%test correct
		go.Experiment.Info.Set('go',{'trial','test','correct'},GenOrder([false true]));
%load some images
	strDirImage	= DirAppend(go.Experiment.File.GetDirectory('code'),'@GridOp','image');
	
	strPathArrow	= PathUnsplit(strDirImage,'arrow','bmp');
	go.arrow		= ind2rgb(uint8(~imread(strPathArrow)),[GO.Param('color','back');GO.Param('color','fore')]);
	
	cOp		= {'cw';'ccw';'h';'v'};
	cPathOp	= cellfun(@(op) PathUnsplit(strDirImage,op,'bmp'),cOp,'uni',false);
	go.op	= cellfun(@(f) ind2rgb(uint8(~imread(f)),[GO.Param('color','back');GO.Param('color','fore')]),cPathOp,'uni',false);
%set the initial reward
	go.reward	= GO.Param('reward','base');

go.Experiment.Info.Set('go','prepared',true);


%------------------------------------------------------------------------------%
function order = GenOrder(k)
	n		= numel(k);
	nRepK	= ceil(nTrialPer/n);
	
	order	= blockdesign(k,nRepK,nRun);
	order	= order(:,1:nTrialPer);
end
%------------------------------------------------------------------------------%

end
