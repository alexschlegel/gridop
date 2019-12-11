% Analysis_20140724_RandomEffects
% address the Science reviewer's concern about our correlation analysis being
% fixed effects but needing to be random effects.  i'm not convinved this is an
% issue, considering we followed the method used by Kriegeskorte and as far as
% i can tell everyone else doing RSA. but, try random-effects analyses to see
% how things go.
nCore	= 1;
strNameAnalysis	= '20140724_randomeffects';

strDirOut	= DirAppend(strDirAnalysis,strNameAnalysis);
CreateDirPath(strDirOut);

ifo			= GO.SubjectInfo;
cSubject	= ifo.code.fmri;
nSubject	= numel(cSubject);

cMask		= GO.UnionMasks;
cMask		= [cMask.core];
nMask		= numel(cMask);
cMaskCore	= [cMask; 'core'];

cMaskLabel		= upper(cMask);
cMaskLabel{5}	= 'LOC';

dimPCAGC		= 10;
dimPCACC		= 50;

colMask	=	[
				0	222	222	%dlpfc
				224	0	224	%fef
				255	0	0	%occ
				224	144	0	%pcu
				64	96	255	%loc
				0	160	0	%ppc
			]/255;

%get the model
	conf	= GO.ConfusionModels;
	conf	= conf(1);

%ROI MVPA
	%load the results of the ROI classification analysis
		[res,stat]	= GO.Analyze.ROIMVPA(...
						'subject'		, cSubject	, ...
						'mask'			, cMaskCore	, ...
						'mindim'		, dimPCAGC	, ...
						'ifo'			, ifo		, ...
						'force_mvpa'	, false		, ...
						'cores'			, nCore		  ...
						);
	
	%recompute the stats with the new jackknife option
		statNew	= MVPAClassifyExtraStats(res,...
					'confusion_model'	, conf	  ...
					);
	
	%remake the confusion correlation plot
		h	= GO.Plot.ConfusionCorrelation(statNew,'outdir',strDirOut,'jackknife',true);

%Cross-classification
	%load the results of the ROI cross-classification analysis
		[res,stat,cMaskPair]	= GO.Analyze.CCMVPA(...
									'subject'	, cSubject	, ...
									'mask'		, cMask		, ...
									'dim'		, dimPCACC	, ...
									'ifo'		, ifo		, ...
									'cores'		, nCore		  ...
									);
		nMaskPair				= size(cMaskPair,1);
	
	%recalculate the stats without the same-same pairs (i just realized i FDR
	%corrected with these in there too, oops.)
	%also add in the jackknife statistics
		cMaskPairDiff	= cMaskPair(1:end-nMask,:);
		nMaskPairDiff	= size(cMaskPairDiff,1);
		
		cScheme	= GO.Param('scheme');
		nScheme	= numel(cScheme);
		
		stat	= struct;
		for kS=1:nScheme
			strScheme	= cScheme{kS};
			
			%accuracy
				stat.(strScheme).acc	= reshape(res.(strScheme).allway.accuracy.mean,[nSubject nMaskPair]);
				stat.(strScheme).acc	= stat.(strScheme).acc(:,1:nMaskPairDiff);
				stat.(strScheme).mAcc	= mean(stat.(strScheme).acc,1)';
				stat.(strScheme).seAcc	= stderr(stat.(strScheme).acc,[],1)';
				
				[h,p,ci,stats]				= ttest(stat.(strScheme).acc,0.25,'tail','right');
				[pThresh,pFDR]				= fdr(p,0.05);
				stat.(strScheme).pAcc		= p';
				stat.(strScheme).pfdrAcc	= pFDR';
				stat.(strScheme).tAcc		= stats.tstat';
				stat.(strScheme).dfAcc		= stats.df';
			
			%confusion
				stat.(strScheme).conf	= reshape(res.(strScheme).allway.confusion,[4 4 nSubject nMaskPair]);
				stat.(strScheme).conf	= stat.(strScheme).conf(:,:,:,1:nMaskPairDiff);
				stat.(strScheme).mConf	= squeeze(mean(stat.(strScheme).conf,3));
				stat.(strScheme).seConf	= squeeze(stderr(stat.(strScheme).conf,[],3));
				
				[r,stats]					= corrcoef2(reshape(conf{1},[],1),reshape(permute(stat.(strScheme).mConf,[3 1 2]),nMaskPairDiff,[]));
				[pThresh,pFDR]				= fdr(stats.p,0.05);
				stat.(strScheme).pConf		= stats.p;
				stat.(strScheme).pfdrConf	= pFDR;
				stat.(strScheme).rConf		= stats.r;
				stat.(strScheme).dfConf		= stats.df;
			
			%jackknifed confusion
				confJK	= permute(stat.(strScheme).conf,[3 1 2 4]);
				for kM=1:nMaskPairDiff
					confCur				= reshape(confJK(:,:,:,kM),nSubject,[]);
					confCurJK			= jackknife(@(x) mean(x,1),confCur);
					confJK(:,:,:,kM)	= reshape(confCurJK,nSubject,4,4);
				end
				stat.(strScheme).confJK		= permute(confJK,[2 3 1 4]);
				
				%these should be the same as non-jackknifed
				stat.(strScheme).mConfJK	= squeeze(mean(stat.(strScheme).confJK,3));
				stat.(strScheme).seConfJK	= squeeze(stderrJK(stat.(strScheme).confJK,[],3));
				
				stat.(strScheme).rConfJK	= NaN(nSubject,nMaskPairDiff);
				for kU=1:nSubject
					for kM=1:nMaskPairDiff
						stat.(strScheme).rConfJK(kU,kM)	= corrcoef2(reshape(conf{1},[],1),reshape(stat.(strScheme).confJK(:,:,kU,kM),1,[]));
					end
				end
				
				stat.(strScheme).mRConfJK	= mean(stat.(strScheme).rConfJK,1)';
				stat.(strScheme).seRConfJK	= stderrJK(stat.(strScheme).rConfJK,[],1)';
				
				[h,p,ci,stats]					= ttestJK(stat.(strScheme).rConfJK,0,0.05,'right',1);
				[pThresh,pFDR]					= fdr(p,0.05);
				stat.(strScheme).pRConfJK		= p';
				stat.(strScheme).pfdrRConfJK	= pFDR';
				stat.(strScheme).tRConfJK		= stats.tstat';
				stat.(strScheme).dfRConfJK		= stats.df';
		end
	
	%what did we get?
		p	= num2cell([stat.shape.pAcc, stat.shape.pfdrAcc, stat.shape.pRConfJK, stat.shape.pfdrRConfJK, stat.operation.pAcc, stat.operation.pfdrAcc, stat.operation.pRConfJK, stat.operation.pfdrRConfJK]);
		lbl	= {'roi1' 'roi2' 'S p/A' 'S pFDR/A' 'S p/C' 'S pFDR/C' 'O p/A' 'O pFDR/A' 'O p/C' 'O pFDR/C'};
		x	= [cMaskPairDiff p]';
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
					
					kPair	= find(strcmp(cMaskPairDiff(:,2),strMaskR) & strcmp(cMaskPairDiff(:,1),strMaskC));
					
					[t(kMR,kMC),t(kMC,kMR)]			= deal(stat.(strScheme).tRConfJK(kPair));
					[p(kMR,kMC),p(kMC,kMR)]			= deal(stat.(strScheme).pRConfJK(kPair));
					[pfdr(kMR,kMC),pfdr(kMC,kMR)]	= deal(stat.(strScheme).pfdrRConfJK(kPair));
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
					'cmax'			, 4.75			, ...
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

%GC MVPA
	%load the results of the Granger-causality classification analysis
		[res,stat] = GO.Analyze.GCMVPA(...
					'subject'	, cSubject	, ...
					'mask'		, cMask		, ...
					'dim'		, dimPCAGC	, ...
					'ifo'		, ifo		, ...
					'cores'		, nCore		  ...
					);
	
	%recompute the stats with the new jackknife option
		stat	= MVPAClassifyExtraStats(res,...
					'confusion_model'	, conf	  ...
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
					
					rJK		= stat.confusion.corr.subjectJK.allway.group.r(kS,kMS,kMD);
					zJK		= fisherz(rJK);
					tJK		= stat.confusion.corr.subjectJK.allway.group.t(kS,kMS,kMD);
					
					p		= stat.confusion.corr.group.allway.p(kS,kMS,kMD);
					pfdr	= stat.confusion.corr.group.allway.pfdr(kS,kMS,kMD);
					
					pJK		= stat.confusion.corr.subjectJK.allway.group.p(kS,kMS,kMD);
					pfdrJK	= stat.confusion.corr.subjectJK.allway.group.pfdr(kS,kMS,kMD);
					
					pacc	= stat.accuracy.p.allway(kS,kMS,kMD);
					pfdracc	= stat.accuracy.pfdr.allway(kS,kMS,kMD);
					
					kSrc	= find(strcmp(strMaskSrc,cMask));
					kDst	= find(strcmp(strMaskDst,cMask));
					
					conn(kSrc,kDst,kS)		= tJK;
					pConn(kSrc,kDst,kS)		= pJK;
					pFDRConn(kSrc,kDst,kS)	= pfdrJK;
					
					%if p <= 0.05 || pacc <= 0.05
					%if pfdr <= 0.05
					if p <= 0.05
						disp(sprintf('  %5s to %5s: t=%.3f p=%.4f/pC=%.4f pjk=%.4f/pjkC=%.4f pacc=%.4f/paccC=%.4f',strMaskSrc,strMaskDst,tJK,p,pfdr,pJK,pfdrJK,pacc,pfdracc));
					elseif pJK <= 0.05
						disp(sprintf('**%5s to %5s: t=%.3f p=%.4f/pC=%.4f pjk=%.4f/pjkC=%.4f pacc=%.4f/paccC=%.4f',strMaskSrc,strMaskDst,tJK,p,pfdr,pJK,pfdrJK,pacc,pfdracc));
					end
				end
			end
		end
		
	%remake the connectivity figures
		cLabel		= cMaskLabel;
		[b,kOrder]	= ismember(cLabel,cMaskLabel);
		colLabel	= colMask(kOrder,:);
		
		col	= GetPlotColors(2);
		col	= col(end:-1:1,:);
		
		%get the line with bounds
			bSize	= pConn<=0.05;
			CSize	= conn(bSize);
			CMin	= min(CSize);
			CMax	= max(CSize);
			
			tLineMin	= 1;
			tLineMax	= 8;
		
		for kS=1:nScheme
			strScheme	= cScheme{kS};
			
			C		= ReorderConfusion(conn(:,:,kS),cMaskLabel,cLabel);
			p		= ReorderConfusion(pConn(:,:,kS),cMaskLabel,cLabel);
			pfdr	= ReorderConfusion(pFDRConn(:,:,kS),cMaskLabel,cLabel);
			
			%bShow	= pfdr<=0.05;
			bShow	= p<=0.05;
			
			%hierarchical network plot
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
					
					e(kE).LineWidth	= MapValue(e(kE).Weight,CMin,CMax,tLineMin,tLineMax);
					e(kE).LineColor	= col(kS,:);
				end
				
				[kSupraR,kSupraC]	= find(p<=0.05 & pfdr>0.05);
				nSupra				= numel(kSupraR);
				disp(sprintf('%s pfdr>0.05:',strScheme));
				for kU=1:nSupra
					kR	= kSupraR(kU);
					kC	= kSupraC(kU);
					disp(sprintf('\t%s to %s: t=%f',cLabel{kR},cLabel{kC},MapValue(C(kR,kC),CMin,CMax,tLineMin,tLineMax)));
				end
				
				g	= biograph.bggui(b);
				f	= get(g.biograph.hgAxes,'Parent');
				
				strPathOut	= PathUnsplit(strDirOut,['gcnetwork-' strScheme],'png');
				saveas(f,PathAddSuffix(strPathOut,'','eps'));
				fig2png(f,strPathOut);
		end
		