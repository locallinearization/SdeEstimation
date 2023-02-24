function FileNames=MakeFileNames(ModEq,ObsEq,ExactSolution)

FileNames{1}=str2func(ModEq);
FileNames{2}=str2func(ObsEq);
if nargin==3
  FileNames{3}=str2func(ExactSolution);
end
