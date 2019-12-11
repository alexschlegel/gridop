function cSubject = Subject(varargin)
% GO.Subject
% 
% Description:	get the subject codes
% 
% Syntax:	cSubject = GO.Subject(<options>)
% 
% In:
% 	<options>:
%		subject:	(<all>) the subject codes to return
%		state:		(<default>) see GO.SubjectInfo
% 
% Out:
% 	cSubject	- the subject codes
% 
% Updated: 2014-03-15
% Copyright 2014 Alex Schlegel (schlegel@gmail.com).  This work is licensed
% under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
% License.
persistent mSubject;

if isempty(mSubject)
	mSubject	= mapping;
end

opt	= ParseArgs(varargin,...
		'subject'	, {}	, ...
		'state'		, []	  ...
		);

if isempty(opt.subject)
	cSubject	= mSubject(opt.state);
	
	if isempty(cSubject)
		s					= GO.SubjectInfo('state',opt.state);
		cSubject			= s.code.fmri;
		mSubject(opt.state)	= cSubject;
	end
else
	cSubject	= ForceCell(opt.subject);
end
