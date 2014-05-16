%FileList = dir('*.nii');
FileList = spm_select(Inf,'image')
NFiles = size(FileList,1);
BaseDir = pwd;
OutDir = BaseDir;
%thr = 1.96;
%k = 20;
bgImage = '/share/studies/CogRes/GroupAnalyses/ModMedCogRes/wmT1_P00004200_S0001.nii';
for i = 1:NFiles
    thrstatImage = FileList(i,:);%fullfile(BaseDir,FileList(i).name);
    % remove the index off the end of the file name
    [PathName FileName Ext] = fileparts(thrstatImage);
    thrstatImage = fullfile(PathName,FileName);
        if ~exist(fullfile(PathName,'ResultFigures'),'dir')
            mkdir(fullfile(PathName,'ResultFigures'));
        end
        RenderstatImage = fullfile(PathName,'Render.img');
        SliceImage = fullfile(PathName,'ResultFigures',[FileName '.png']);
        %if exist(statImage)
        threshold stat image
        str = ['fslmaths ' statImage ' -thr ' num2str(thr) ' '  thrstatImage];
        [a b] = unix(str);
        str = ['overlay 0 0 ' bgImage ' -a ' thrstatImage ' ' num2str(thr) ' 7 ' RenderstatImage ' y']
        unix(str);
        str = ['fslstats ' RenderstatImage ' -R']
        [a b] = unix(str);
        fSpace = findstr(b,' ');
        lower = b(1:fSpace(1)-1);
        upper = b(fSpace(1)+1:end-1);
        %str = ['slicer ' RenderstatImage ' -l render1 -i ' lower ' 7 -z -35 ' SliceImage]
        str = ['slicer ' RenderstatImage ' -l render1 -i ' lower ' 7 -a ' SliceImage ]
        unix(str);
        

        % Create POS  versions
        [PathName FileName Ext] = fileparts(statImage);
        SliceImage = fullfile(PathName,'ResultFigures',['POS_' FileName '_T' num2str(thr) '_k' num2str(k) '.png']);
        RenderstatImage = fullfile(PathName,'Render.img');
                %if exist(statImage)
        str = ['overlay 0 0 ' bgImage ' -a ' thrstatImage ' ' num2str(thr) ' 7 ' RenderstatImage ' y']
        unix(str);
        str = ['fslstats ' RenderstatImage ' -R']
        [a b] = unix(str);
        fSpace = findstr(b,' ');
        lower = b(1:fSpace(1)-1);
        upper = b(fSpace(1)+1:end-1);
        %str = ['slicer ' RenderstatImage ' -l render1 -i ' lower ' 7 -z -35 ' SliceImage]
        str = ['slicer ' RenderstatImage ' -l render1 -i ' lower ' 7 -a ' SliceImage ]
        unix(str);
        
        % Create NEG  versions
        [PathName FileName Ext] = fileparts(statImage);
        
        SliceImage = fullfile(PathName,'ResultFigures',['NEG_' FileName '_T' num2str(thr) '_k' num2str(k) '.png']);
        thrstatImage = fullfile(PathName,'tempStat.nii');
        RenderstatImage = fullfile(PathName,'Render.img');
                %if exist(statImage)
        % threshold stat image
        str = ['fslmaths ' statImage ' -mul -1 -thr ' num2str(thr) ' '  thrstatImage];
        [a b] = unix(str);
        str = ['overlay 0 0 ' bgImage ' -a ' thrstatImage ' ' num2str(thr) ' 7 ' RenderstatImage ' y']
        unix(str);
        str = ['fslstats ' RenderstatImage ' -R']
        [a b] = unix(str);
        fSpace = findstr(b,' ');
        lower = b(1:fSpace(1)-1);
        upper = b(fSpace(1)+1:end-1);
        %str = ['slicer ' RenderstatImage ' -l render1 -i ' lower ' 7 -z -35 ' SliceImage]
        str = ['slicer ' RenderstatImage ' -l rendersea -i ' lower ' 7 -a ' SliceImage ]
        unix(str);
    end
    
end
