function n = NumTrials(ifo)
% GO.Paper.NumTrials
% 
% Description:	calculate the number of correct trials for each condition
% 
% Syntax:	n = GO.Paper.NumTrials(ifo)
% 
% Updated: 2014-07-18
% Copyright 2014 Alex Schlegel (schlegel@gmail.com).  This work is licensed
% under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
% License.
nSubject	= size(ifo.correct,1);

cScheme	= {'operation','shape','shapeop'};
nScheme	= numel(cScheme);

for kS=1:nScheme
	strScheme	= cScheme{kS};
	
	kCondition	= unique(ifo.(strScheme));
	nCondition	= numel(kCondition);
	
	n.(strScheme)	= NaN(nSubject,nCondition);
	
	correct	= reshape(ifo.correct,nSubject,[]);
	trial	= reshape(ifo.(strScheme),nSubject,[]);
	for kC=1:nCondition
		n.(strScheme)(:,kC)	= sum(correct .* trial==kCondition(kC),2);
	end
	
	n.stat.(strScheme).all.m	= mean(n.(strScheme)(:),1);
	n.stat.(strScheme).all.se	= stderr(n.(strScheme)(:),[],1);
	n.stat.(strScheme).m		= mean(n.(strScheme),1);
	n.stat.(strScheme).se		= stderr(n.(strScheme),[],1);
end

n.stat.all.m	= mean(reshape([n.shape n.operation],[],1));
n.stat.all.se	= stderr(reshape([n.shape n.operation],[],1));
