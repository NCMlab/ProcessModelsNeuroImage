function  [xSPM hReg] = SPMDisplayThresholdImage(Input, HeightThr, ExtentThr)

% Filename = SPMDisplayThresholdImage.m
% Dependent on:  DisplayCov.m
%                 jrs_results_ui.m
% Written by: Jason Steffener
% Date: 09/07
%

clear mni tal value cluster;


  
warning('off','all')
    SCCSid = '$Rev: 816 $';
    
    
   %-Initialise
    %----------------------------------------------------------------------
    SPMid      = spm('FnBanner',mfilename,SCCSid);
    [Finter,Fgraph,CmdLine] = spm('FnUIsetup','Stats: Results');
    FS         = spm('FontSizes');
    
    % clear satfig if it exists
    %----------------------------------------------------------------------
    hSat       = findobj('tag','Satellite');
    spm_figure('clear',hSat);

    %-Get thresholded xSPM data and parameters of design
    %=======================================================================
    if nargin == 3
        [SPM xSPM] = DisplayCov(Input, HeightThr, ExtentThr);
    elseif nargin == 1
        [SPM xSPM] = DisplayCov(Input);
    else
        [SPM xSPM] = DisplayCov;
    end
%	[SPM,xSPM] = spm_getSPM;


    if isempty(xSPM) 
	varargout = {[],[],[]};
	return;
    end

    M         = SPM.xVol.M;
    DIM       = SPM.xVol.DIM;
    try
        units = SPM.xVol.units;
    catch
        units = {'mm' 'mm' 'mm'};
    end
    % ensure pwd = swd so that relative filenames are valid
    %----------------------------------------------------------------------
    
    cd(SPM.swd)

    %-Setup Results User Interface; Display MIP, design matrix & parameters
    %======================================================================
    spm('FigName',['SPM{',xSPM.STAT,'}: Results'],Finter,CmdLine);




    %-Setup results GUI
    %----------------------------------------------------------------------
    spm_figure('Clear',Finter)
    hReg      = jrs_results_ui('SetupGUI',M,DIM,xSPM,Finter);

    %-Setup design interrogation menu
    %----------------------------------------------------------------------
%     hDesRepUI = spm_DesRep('DesRepUI',SPM);
%     figure(Finter)

    %-Setup Maximium intensity projection (MIP) & register
    %----------------------------------------------------------------------
    %spm('defaults','FMRI')
    hMIPax = axes('Parent',Fgraph,'Position',[0.05 0.60 0.55 0.36],'Visible','off');
    hMIPax = spm_mip_ui(xSPM.Z,xSPM.XYZmm,M,DIM,hMIPax,units);

    spm_XYZreg('XReg',hReg,hMIPax,'spm_mip_ui');
    if xSPM.STAT == 'P'
        str = xSPM.STATstr;
    else
        str = ['SPM\{',xSPM.STATstr,'\}'];
    end
    text(240,260,str,...
        'Interpreter','TeX',...
        'FontSize',FS(14),'Fontweight','Bold',...
        'Parent',hMIPax)


    %-Print comparison title
    %----------------------------------------------------------------------
    hTitAx = axes('Parent',Fgraph,...
        'Position',[0.02 0.95 0.96 0.02],...
        'Visible','off');

    text(0.5,0,xSPM.title,'Parent',hTitAx,...
        'HorizontalAlignment','center',...
        'VerticalAlignment','baseline',...
        'FontWeight','Bold','FontSize',FS(14))


    %-Print SPMresults: Results directory & thresholding info
    %----------------------------------------------------------------------
    hResAx = axes('Parent',Fgraph,...
        'Position',[0.05 0.55 0.45 0.05],...
        'DefaultTextVerticalAlignment','baseline',...
        'DefaultTextFontSize',FS(9),...
        'DefaultTextColor',[1,1,1]*.7,...
        'Units','points',...
        'Visible','off');
    AxPos = get(hResAx,'Position'); set(hResAx,'YLim',[0,AxPos(4)])
    h     = text(0,24,'SPMresults:','Parent',hResAx,...
        'FontWeight','Bold','FontSize',FS(14));
    text(get(h,'Extent')*[0;0;1;0],24,spm_str_manip(SPM.swd,'a30'),'Parent',hResAx)
    try
        thresDesc = xSPM.thresDesc;
        text(0,12,sprintf('Height threshold %c = %0.6f  {%s}',xSPM.STAT,xSPM.u,thresDesc),'Parent',hResAx)
    catch
        text(0,12,sprintf('Height threshold %c = %0.6f',xSPM.STAT,xSPM.u),'Parent',hResAx)
    end
    text(0,00,sprintf('Extent threshold k = %0.0f voxels',xSPM.k), 'Parent',hResAx)

    %-Store handles of results section Graphics window objects
    %----------------------------------------------------------------------
    H  = get(Fgraph,'Children');
    H  = findobj(H,'flat','HandleVisibility','on');
    H  = findobj(H);
    Hv = get(H,'Visible');
    set(hResAx,'Tag','PermRes','UserData',struct('H',H,'Hv',{Hv}))

    %-Finished results setup
    %----------------------------------------------------------------------
    varargout = {hReg,xSPM,SPM};
    spm('Pointer','Arrow')
    % spm_list in SPM8 needs this field for FDR values
    if ~isfield(xSPM,'uc')
        xSPM.uc = [];
    end
    
    

   
    if spm('ver') == 'SPM8'
    % There is an issue with calculation of the FDR values with SPM8.
    % To circomnavigate this the defaults need to be changed.
        OrigDefault = spm_get_defaults('stats.topoFDR');
        spm_get_defaults('stats.topoFDR',0);
        spm_get_defaults('stats.topoFDR')
    end
    
     spm_list('list',xSPM,hReg)
    if spm('ver') == 'SPM8'
    % reset the default value
        spm_get_defaults('stats.topoFDR',OrigDefault);
    end
        
%%% Little add-on by Chris

clear *mni tal *cluster *value index;
temp=ans;





% for i=1:size(temp.dat,1)
%   mmni(i,:)=temp.dat{i,end};   %MNI coordinates
%   mvalue(i)=temp.dat{i,end-3};     %voxel values
% 
%   if isempty(temp.dat{i,4})
%      temp.dat{i,4}=-1;
%   end
%   
%   mcluster(i)=temp.dat{i,4};
% end
% 
% 
% %threshold with cluster value chosen
% index=find(mcluster>=xSPM.k);
% 
% 
% mni=mmni(index,:);
% 
% 
% tal=round(mni2tal(mni))
% value=mvalue(index)
% cluster=mcluster(index)
% 


