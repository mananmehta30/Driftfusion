% designate subfolber to save to
subfolder = 'B8_Ag-MAPI-Al_mobcat\sc_l-sc_r';


% base folder is Notes folder on OneDrive
FolderDir = strcat('C:\Users\felix\OneDrive\Dokumente\_College\Project\Notes\'...
    ,subfolder);

% Create subfolder if it doesn't exist
if ~exist(FolderDir,'dir')
    disp('Subfolder doesn''t exist, creating folder')
    mkdir(FolderDir)
end


FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
%%
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = num2str(get(FigHandle, 'Number'));
  set(0, 'CurrentFigure', FigHandle);
  
  if exist([FolderDir,'\',FigName,'.png'],'file') || exist([FolderDir,'\',FigName,'.fig'],'file')
      break
  else
    savefig(fullfile(FolderDir, [FigName '.fig']));
    saveas(FigHandle,fullfile(FolderDir, [FigName '.png']));
  end
end