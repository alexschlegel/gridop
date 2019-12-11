function s = UnionMasks()
% GO.UnionMasks
% 
% Description:	get a struct of masks to include in each type of union
% 
% Syntax:	s = GO.UnionMasks()
% 
% Updated: 2014-03-03
% Copyright 2014 Alex Schlegel (schlegel@gmail.com).  This work is licensed
% under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
% License.

s			= struct;
s.network	= {'cere','dlpfc','fef','fo','mfc','mtl','occ','pcu','pitc','ppc','sef','thal'}';
s.core		= {'dlpfc','fef','occ','pcu','pitc','ppc'}';
s.corethal	= [s.core; 'thal'];

cBase	= fieldnames(s);
nBase	= numel(cBase);

for kB=1:nBase
	strBase	= cBase{kB};
	
	cMask	= s.(strBase);
	nMask	= numel(cMask);
	
	for kM=1:nMask
		strMask	= cMask{kM};
		
		s.([strBase '_no_' strMask])	= setdiff(cMask,strMask); 
	end
end
