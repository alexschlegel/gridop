function [n,cSubject] = TrialsToCriterion(varargin)
% GO.TrialsToCriterion
% 
% Description:	calculate the number of trials to criterion during the practice
%				session for each subject
% 
% Syntax:	[n,cSubject] = TrialsToCriterion(<options>)
% 
% In:
% 	<options>:
%		ifo:	(<load>) the subject info struct (GO.SubjectInfo)
% 
% Out:
% 	n			- number of trials until criterion was reached for each subject
%	cSubject	- the subject ids
% 
% Updated: 2014-05-29
% Copyright 2014 Alex Schlegel (schlegel@gmail.com).  This work is licensed
% under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
% License.
opt	= ParseArgs(varargin,...
		'ifo'	, []	  ...
		);

if isempty(opt.ifo)
	opt.ifo	= GO.SubjectInfo;
end

cSubject	= opt.ifo.id;
nSubject	= numel(cSubject);

n	= NaN(nSubject,1);

for kS=1:nSubject
	bHistory	= opt.ifo.practice.history{kS};
	nHistory	= numel(bHistory);
	
	for kH=10:nHistory
		if all(bHistory(kH-10 + (1:10)))
			n(kS)	= kH;
			break;
		end
	end
end
