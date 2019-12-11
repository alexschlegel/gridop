function Behavioral(varargin)
% GO.Analyze.Behavioral
% 
% Description:	behavioral analyses
% 
% Syntax:	GO.Analyze.Behavioral()
% 
% Updated: 2014-05-21
% Copyright 2014 Alex Schlegel (schlegel@gmail.com).  This work is licensed
% under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
% License.
ifo	= GO.SubjectInfo;

%mean performance accuracy
	m	= 100*mean(ifo.correct(:));
	
	disp(sprintf('mean accuracy: %.2f%%',m));

%difficulty comparisons
	[p,anovatab,stats]	= anova1(ifo.mean.shape.correct,[],'off');