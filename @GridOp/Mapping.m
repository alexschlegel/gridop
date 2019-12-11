function Mapping(go,varargin)
% GridOp.Mapping
% 
% Description:	show the stimulus and operation mappings
% 
% Syntax:	go.Mapping(<options>)
%
% In:
%	<options>:
%		wait:	(true) true to wait for user input before returning
% 
% Updated: 2013-09-25
% Copyright 2013 Alex Schlegel (schlegel@gmail.com).  This work is licensed
% under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
% License.
opt	= ParseArgs(varargin,...
		'wait'	, true	  ...
		);

%get the stimulus and operation mappings for the current subject
	mapStim	= go.Experiment.Subject.Get('map_stim');
	mapOp	= unless(go.Experiment.Subject.Get('map_op'),[1;3;2;4]);
%get the stimuli
	cStim	= arrayfun(@(k) GO.Stim.Stimulus(k,'map',mapStim),(1:4)','uni',false);

%open a texture
	sTexture	= switch2(go.Experiment.Info.Get('go','session'),1,800,2,1000);
	go.Experiment.Window.OpenTexture('mapping',[sTexture sTexture]);
%show the stimuli
	go.Experiment.Show.Text('<size:1><color:marigold>shapes</color></size>',[0 -3.5],'window','mapping');
	
	strStim	= GO.Param('prompt','stimulus');
	
	for k=1:4
		go.Experiment.Show.Image(cStim{k},[4*(k-1)-6 -1.75],2.5,'window','mapping');
		go.Experiment.Show.Text(['<size:1><style:normal>' strStim(k) '</style></size>'],[4*(k-1)-6 0.5],'window','mapping');
	end

%show the operations
	go.Experiment.Show.Text('<size:1><color:marigold>operations</color></size>',[0 2.25],'window','mapping');
	
	strOp	= GO.Param('prompt','operation');
	
	for k=1:4
		go.Experiment.Show.Image(go.op{mapOp(k)},[4*(k-1)-6 4],2.5,'window','mapping');
		go.Experiment.Show.Text(['<size:1><style:normal>' strOp(k) '</style></size>'],[4*(k-1)-6 6.25],'window','mapping');
	end

%show the instructions screen
	if opt.wait
		fResponse	= [];
		strPrompt	= [];
	else
		fResponse	= false;
		strPrompt	= ' ';
	end
	
	fResponse	= conditional(opt.wait,[],false);
	go.Experiment.Show.Instructions('',...
					'figure'	, 'mapping'	, ...
					'fresponse'	, fResponse	, ...
					'prompt'	, strPrompt	  ...
					);

%remove the texture
	go.Experiment.Window.CloseTexture('mapping');
	
