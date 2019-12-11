% Analysis_20140428_All
% Compilation of all final analyses:
% 1) ROI MVPA
%	a. pre-whiten data using PCA
%	b. classify shape and operation within this component space
%	c. accuracies + correlate confusion matrices with model
% 2) GC MVPA
%	a. take the top 10 PCA components for each ROI and calculate GC patterns
%	   between ROI pairs
%	b. classify shape + operation using GC patterns for each pair and direction
%	c. accuracies + correlate confusion matrices with model
% 3) ROI cross-classification
%	a. take the top 50 PCA components for each ROI.
%	b. for each pair of ROIs, match up these 50 components so that the ROIs are
%	   now ideally in the same space.
%	c. cross-classify between the masks
nCore			= 12;
strNameAnalysis	= '20140428_all';

strDirOut	= DirAppend(strDirAnalysis,strNameAnalysis);
CreateDirPath(strDirOut);

ifo			= GO.SubjectInfo;
cSubject	= ifo.code.fmri;
cMask		= GO.UnionMasks;
cMask		= [cMask.core];
nMask		= numel(cMask);
cMaskCore	= [cMask; 'core'];

cMaskLabel		= upper(cMask);
cMaskLabel{5}	= 'LOC';

dimPCAGC	= 10;
dimPCACC	= 50;

colMask	=	[
				0	222	222	%dlpfc
				224	0	224	%fef
				255	0	0	%occ
				224	144	0	%pcu
				64	96	255	%loc
				0	160	0	%ppc
			]/255;

%open a MATLAB pool
	MATLABPoolOpen(nCore);

%ROI Classification!
	[res,stat]	= GO.Analyze.ROIMVPA(...
					'subject'		, cSubject	, ...
					'mask'			, cMaskCore	, ...
					'mindim'		, dimPCAGC	, ...
					'ifo'			, ifo		, ...
					'cores'			, nCore		  ...
					);

	%figures
		h	= GO.Plot.Accuracy(stat,'outdir',strDirOut);
		h	= GO.Plot.ConfusionCorrelation(stat,'outdir',strDirOut);
		h	= GO.Plot.Confusion(stat,'outdir',strDirOut);

%GC Connectivity Pattern Classification!
	[res,stat] = GO.Analyze.GCMVPA(...
					'subject'	, cSubject	, ...
					'mask'		, cMask		, ...
					'dim'		, dimPCAGC	, ...
					'ifo'		, ifo		, ...
					'cores'		, nCore		  ...
					);
	
	%what did we get?
		cScheme	= stat.label{1};
		nScheme	= numel(cScheme);
		
		[conn,pConn,pFDRConn]	= deal(NaN(nMask,nMask,nScheme));
		
		for kS=1:nScheme
			strScheme	= stat.label{1}{kS};
			
			disp(sprintf('%s:',strScheme));
			
			cMaskSrc	= stat.label{2};
			nMaskSrc	= numel(cMaskSrc);
			
			for kMS=1:nMaskSrc
				strMaskSrc	= cMaskSrc{kMS};
				
				cMaskDst		= cMaskSrc;
				bSelf			= strcmp(strMaskSrc,cMaskDst);
				cMaskDst(bSelf)	= [];
				nMaskDst		= numel(cMaskDst);
				
				for kMD=1:nMaskDst
					strMaskDst	= cMaskDst{kMD};
					
					r		= stat.confusion.corr.group.allway.r(kS,kMS,kMD);
					z		= fisherz(r);
					
					p		= stat.confusion.corr.group.allway.p(kS,kMS,kMD);
					pfdr	= stat.confusion.corr.group.allway.pfdr(kS,kMS,kMD);
					
					pacc	= stat.accuracy.p.allway(kS,kMS,kMD);
					paccfdr	= stat.accuracy.pfdr.allway(kS,kMS,kMD);
					
					kSrc	= find(strcmp(strMaskSrc,cMask));
					kDst	= find(strcmp(strMaskDst,cMask));
					
					conn(kSrc,kDst,kS)		= z;
					pConn(kSrc,kDst,kS)		= p;
					pFDRConn(kSrc,kDst,kS)	= pfdr;
					
					%if p <= 0.05 || pacc <= 0.05
					%if pfdr <= 0.05
					if p <= 0.05
						disp(sprintf('   %6s to %6s: r=%f, p=%f, pfdr=%f, pacc=%f, paccfdr=%f',strMaskSrc,strMaskDst,r,p,pfdr,pacc,paccfdr));
					end
				end
			end
		end
	
	%connectivity figures
		%custom anterior to posterior order (nope, i guess i decided against this)
			cLabel		= cMaskLabel;
			[b,kOrder]	= ismember(cLabel,cMaskLabel);
			colLabel	= colMask(kOrder,:);
		
		col	= GetPlotColors(2);
		col	= col(end:-1:1,:);
		
		for kS=1:nScheme
			strScheme	= cScheme{kS};
			
			C		= ReorderConfusion(conn(:,:,kS),cMaskLabel,cLabel);
			p		= ReorderConfusion(pConn(:,:,kS),cMaskLabel,cLabel);
			pfdr	= ReorderConfusion(pFDRConn(:,:,kS),cMaskLabel,cLabel);
			
			%hierarchical network plot
				%get the line with bounds
					rAll	= stat.confusion.corr.group.allway.r;
					zAll	= fisherz(rAll);
					pAll	= stat.confusion.corr.group.allway.p;
					pfdrAll	= stat.confusion.corr.group.allway.pfdr;
					
					bShow	= pfdr<=0.05;
					bSize	= pAll<=0.05;
					
					zSize	= zAll(bSize);
					zMin	= min(zSize);
					zMax	= max(zSize);
					
					tLineMin	= 1;
					tLineMax	= 8;
				
				x			= C;
				x(~bShow)	= 0;
				x(isnan(x))	= 0;
				
				b		= biograph(x,cLabel);
				
				n		= b.Nodes;
				nNode	= numel(n);
				for kN=1:nNode
					n(kN).Color		= colMask(kN,:);
					n(kN).LineColor	= colMask(kN,:);
					n(kN).TextColor	= GetGoodTextColor(n(kN).Color);
				end
				
				e		= b.Edges;
				nEdge	= numel(e);
				for kE=1:nEdge
					sNode	= regexp(e(kE).ID,'(?<src>.+) -> (?<dst>.+)','names');
					kSrc	= find(strcmp(cMask,sNode.src));
					kDst	= find(strcmp(cMask,sNode.dst));
					
					e(kE).LineWidth	= MapValue(e(kE).Weight,zMin,zMax,tLineMin,tLineMax);
					e(kE).LineColor	= col(kS,:);
				end
				
				[kSupraR,kSupraC]	= find(p<=0.05 & pfdr>0.05);
				nSupra				= numel(kSupraR);
				disp(sprintf('%s pfdr>0.05:',strScheme));
				for kU=1:nSupra
					kR	= kSupraR(kU);
					kC	= kSupraC(kU);
					disp(sprintf('\t%s to %s: t=%f',cLabel{kR},cLabel{kC},MapValue(C(kR,kC),zMin,zMax,tLineMin,tLineMax)));
				end
				
				g	= biograph.bggui(b);
				f	= get(g.biograph.hgAxes,'Parent');
				
				strPathOut	= PathUnsplit(strDirOut,['gcnetwork-' strScheme],'png');
				saveas(f,PathAddSuffix(strPathOut,'','eps'));
				fig2png(f,strPathOut);
		end
		
%Cross Classification!
	[res,stat,cMaskPair] = GO.Analyze.CCMVPA(...
							'subject'	, cSubject	, ...
							'mask'		, cMask		, ...
							'dim'		, dimPCACC	, ...
							'ifo'		, ifo		, ...
							'cores'		, nCore		  ...
							);
	
	%what did we get?
		kModel	= 1;
		
		p	= num2cell([stat.shape.pAcc, stat.shape.pfdrAcc, stat.shape.modelcompare(kModel).pConf, stat.shape.modelcompare(kModel).pfdrConf, stat.operation.pAcc, stat.operation.pfdrAcc, stat.operation.modelcompare(kModel).pConf, stat.operation.modelcompare(kModel).pfdrConf]);
		lbl	= {'roi1' 'roi2' 'S p/A' 'S pFDR/A' 'S p/C' 'S pFDR/C' 'O p/A' 'O pFDR/A' 'O p/C' 'O pFDR/C'};
		x	= [cMaskPair p]';
		disp(sprintf(['%5s\t%5s\t' repmat('%8s\t',[1 size(p,2)])],lbl{:}));
		disp(sprintf(['%5s\t%5s\t' repmat('%.6f\t',[1 size(p,2)]) '\n'],x{:}));
	
	%figures!
		%custom order
			cLabel	=		upper({
								'dlpfc'
								'ppc'
								'occ'
								'loc'
								'pcu'
								'fef'
							});
		
		cScheme	= fieldnames(stat);
		nScheme	= numel(cScheme);
		
		col	= GetPlotColors(2);
		
		for kS=1:nScheme
			strScheme	= cScheme{kS};
			
			[t,p,pfdr]	= deal(NaN(nMask));
			
			for kMR=1:nMask
				strMaskR	= cMask{kMR};
				
				for kMC=1:kMR-1
					strMaskC	= cMask{kMC};
					
					kPair	= find(strcmp(cMaskPair(:,2),strMaskR) & strcmp(cMaskPair(:,1),strMaskC));
					
					[t(kMR,kMC),t(kMC,kMR)]			= deal(stat.(strScheme).tAcc(kPair));
					[p(kMR,kMC),p(kMC,kMR)]			= deal(stat.(strScheme).pAcc(kPair));
					[pfdr(kMR,kMC),pfdr(kMC,kMR)]	= deal(stat.(strScheme).pfdrAcc(kPair));
				end
			end
			
			t		= ReorderConfusion(t,cMaskLabel,cLabel);
			p		= ReorderConfusion(p,cMaskLabel,cLabel);
			pfdr	= ReorderConfusion(pfdr,cMaskLabel,cLabel);
			
			h	= alexplot(t,...
					'sig'			, p				, ...
					'sigcorr'		, pfdr			, ...
					'label'			, cLabel		, ...
					'lut'			, col(kS,:)		, ...
					'colorbar'		, false			, ...
					'ring_radius'	, 0.8			, ...
					'ring_phase'	, pi/6			, ...
					'cmin'			, 1.75			, ...
					'cmax'			, 4				, ...
					'arcmethod'		, 'line'		, ...
					'arcwidth'		, 'scale'		, ...
					'arcwidthmin'	, 1				, ...
					'arcwidthmax'	, 12			, ...
					'type'			, 'connection'	  ...
					);
			
			strPathOut	= PathUnsplit(strDirOut,['crossclassification-' strScheme],'png');
			fig2png(h.hF,strPathOut);
		end
		
		%to get line widths
			h	= alexplot([NaN NaN NaN; 1 NaN NaN; 2 2 NaN],...
					'colorbar'		, false			, ...
					'arcmethod'		, 'line'		, ...
					'arcwidth'		, 'scale'		, ...
					'arcwidthmin'	, 1				, ...
					'arcwidthmax'	, 12			, ...
					'type'			, 'connection'	  ...
					);
			
			strPathOut	= PathUnsplit(strDirOut,'crossclassification-widths','png');
			fig2png(h.hF,strPathOut);

%close the MATLAB pool
	MATLABPoolClose;
