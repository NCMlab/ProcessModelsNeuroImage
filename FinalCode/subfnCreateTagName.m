function Tag = subfnCreateTagName(AllData)

switch AllData.ModelNum
    case '1'
        Tag = sprintf('Model%s_X%s_M%s_Y%s_COV',AllData.ModelNum,AllData.names.X,AllData.names.M{1},AllData.names.Y);
        for i = 1:length( AllData.names.COV)
            Tag = sprintf('%s%s_',Tag,AllData.names.COV{i});
        end
        Tag = sprintf('%sNboot%d',Tag,AllData.Nboot);

    case '4'
        Tag = sprintf('Model%s_X%s_M%s_Y%s_COV',AllData.ModelNum,AllData.names.X,AllData.names.M{1},AllData.names.Y);
        for i = 1:length( AllData.names.COV)
            Tag = sprintf('%s%s_',Tag,AllData.names.COV{i});
        end
        Tag = sprintf('%sNboot%d',Tag,AllData.Nboot);
    case '6'
        Tag = sprintf('Model%s_X%s_M%s_M%s_Y%s_COV',AllData.ModelNum,AllData.names.X,AllData.names.M{1},AllData.names{2},AllData.names.Y);
        for i = 1:length( AllData.names.COV)
            Tag = sprintf('%s%s_',Tag,AllData.names.COV{i});
        end
        Tag = sprintf('%sNboot%d',Tag,AllData.Nboot);    
    case '7'
        Tag = sprintf('Model%s_X%s_M%s_Y%s_W%s_COV',AllData.ModelNum,AllData.names.X,AllData.names.M{1},AllData.names.Y,AllData.names.W);
        for i = 1:length( AllData.names.COV)
            Tag = sprintf('%s%s_',Tag,AllData.names.COV{i});
        end
        Tag = sprintf('%sNboot%d',Tag,AllData.Nboot);
    case '14'
        Tag = sprintf('Model%s_X%s_M%s_Y%s_V%s_COV',AllData.ModelNum,AllData.names.X,AllData.names.M{1},AllData.names.Y,AllData.names.V);
        for i = 1:length( AllData.names.COV)
            Tag = sprintf('%s%s_',Tag,AllData.names.COV{i});
        end
        Tag = sprintf('%sNboot%d',Tag,AllData.Nboot);
    case '58'
        Tag = sprintf('Model%s_X%s_M%s_Y%s_W%s_V%s_COV',AllData.ModelNum,AllData.names.X,AllData.names.M{1},AllData.names.Y,AllData.names.W,AllData.names.V);
        for i = 1:length( AllData.names.COV)
            Tag = sprintf('%s%s_',Tag,AllData.names.COV{i});
        end
        Tag = sprintf('%sNboot%d',Tag,AllData.Nboot);
    case '74'
        Tag = sprintf('Model%s_X%s_M%s_Y%s_COV',AllData.ModelNum,AllData.names.X,AllData.names.M{1},AllData.names.Y);
        for i = 1:length( AllData.names.COV)
            Tag = sprintf('%s%s_',Tag,AllData.names.COV{i});
        end
        Tag = sprintf('%sNboot%d',Tag,AllData.Nboot);
    case '75'
        Tag = sprintf('Model%s_X%s_M%s_M%s_Y%s_COV',AllData.ModelNum,AllData.names.X,AllData.names.M{1},AllData.names{2},AllData.names.Y);
        for i = 1:length( AllData.names.COV)
            Tag = sprintf('%s%s_',Tag,AllData.names.COV{i});
        end
        Tag = sprintf('%sNboot%d',Tag,AllData.Nboot);
end