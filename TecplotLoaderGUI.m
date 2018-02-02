function TecplotLoaderGUI
clc
clear
% TECPLOTLOADERGUI Open and plot Tecplot-like files
% Bu
% Bu
% Julio Barros - 2015

% Create the UI
ui = figure('Visible','off','Position',[10,10,1024,780]);

loadFile = uicontrol('Style','pushbutton',...
    'String','Load Vector File','Position',[315,220,70,25],...
    'Callback',@loadFile_Callback);

ha = axes('Units','Pixels','Position',[40,60,640,480]);
align(loadFile,'Center','None');

% Initialize the UI.
ui.Units = 'normalized';
ha.Units = 'normalized';
loadFile.Units = 'normalized';

% Move the window to the center of the screen.
movegui(ui,'center')

%Make the UI visible.
ui.Visible = 'on';

    function velData = loadFile_Callback(source,data)
        if ispc == 1
            slash = '\';
        else
            slash = '/';
        end
        
        [fileName,pathName] = uigetfile({'*.dat';'*.vec';'*.*'},'Select File');
        
        % Open Tecplot file and load
        fname = [pathName slash fileName];
        velData = importdata(fname);
        velData.source = fname;
        
        header = velData.textdata{:};
        res = regexp(header,'VARIABLES *= *(?<vars>.*) *', 'names');
        if ~isempty(res),
            res = regexp([res(1).vars ','],'[" ]*(?<name>[^,]*?)[" ]*,', 'tokens');
            for ii = 1:length(res),
                
                %oname{ii} = res{ii}{1};
                vname{ii} = regexprep(res{ii}{1},'[W]+','_');
                
                if ~isempty(strfind(vname{ii},'CHC'))
                    break
                end
                
            end
            velData.vars = vname;
        end
        N = length(vname);
        
        % Read Zone header
        res = regexp(header,'I=(?<I>[^"]+), J=(?<J>[^"]+), K=(?<K>[^"]+), F=(?<F>[^"]+)', 'names');
        velData.I = str2double(res.I);
        velData.J = str2double(res.J);
        velData.K = str2double(res.K);
        
    end

end