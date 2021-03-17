function collectcode(startpath,recurse,delefirst)
%
% COLLECTCODE(STARTPATH,RECURSE,DELEFIRST) locates all m-files within the startpath (recursively,
% or non-recursively), and creates a text file containing the complete text of all located m-files.
% (This file is named 'codetext.txt', and is written to the startpath.)
% 
% When prompted, the program will delete all existing versions of 'codetext.txt' on the search path
% before writing a new file of that name.
% 
% The function calls the embedded subfunction 'filesearch' to locate the files.
%
% Requires no input arguments; queries user via uicontrols for necessary information. However,
% valid options include 0,1,2, or 3 arguments, with startpath, recurse (= 0 or 1), and delefirst
% (= 0 or 1) as inputs.
%
% EXAMPLES:
%
% COLLECTCODE prompts the user via uicontrols for all options, including:
%             starting path, inclusion or exclusion of subdirectories, deletion of existing output file
% COLLECTCODE('C:\Mfiles') searches specified directory for m-files; prompts the user for:
%             inclusion or exclusion of subdirectories, deletion of existing output file
% COLLECTCODE('C:\Mfiles',1) searches specified directory recursively for m-files; prompts the user for:
%             deletion of existing output file
% COLLECTCODE('C:\Mfiles',0,1) searches specified directory non-recursively for m-files; deletes existing
%             instances of 'codetext.txt' in the search path.

% Copyright Brett Shoelson, Ph.D.
% Shoelson Consulting
% brett.shoelson@joslin.harvard.edu
% V1: 3/4/99.
% V2: 5/30/01. Extensive modifications include on/off of recursion, handling of different numbers of input arguments.

if ~nargin
	% Get starting path
	[filename,startpath]=uigetfile('*.m','Select any m-file in the desired starting pathname.');
	%Remove trailing '\' character
	if strcmp(startpath(end),'\'),startpath=startpath(1:end-1);end
	% Query for inclusion of subdirectories
	recurse=questdlg('Include subdirectories?','','YES','No','YES');
	if strcmp(recurse,'YES'),recurse=1;else;recurse=0;end
	% Query for deletion of existing codetext.txt file
	delefirst=questdlg('Clear _codetext_ first?','Collectcode.m','YES','No','YES');
	if strcmp(delefirst,'YES'), delefirst=1;else;delefirst=0;end
elseif nargin == 1
	recurse=questdlg('Include subdirectories?','','YES','No','YES');
	if strcmp(recurse,'YES'),recurse=1;else;recurse=0;end
	delefirst=questdlg('Clear _codetext_ first?','Collectcode.m','YES','No','YES');
	if strcmp(delefirst,'YES'), delefirst=1;else;delefirst=0;end
elseif nargin == 2
	delefirst=questdlg('Clear _codetext_ first?','Collectcode.m','YES','No','YES');
	if strcmp(delefirst,'YES'), delefirst=1;else;delefirst=0;end
end

% Delete instances of 'codetext.txt' in search path.
if delefirst
	files=filesearch('txt',startpath,recurse);
	for i=1:length(files)
		if ~isempty(findstr(files{i},'codetext.txt'))
			delete(files{i});
		end
	end
end

% Compile filenames (calls subfunction filesearch)
files=filesearch('m',startpath,recurse);
cd(startpath);

% Open codetext.txt for appending. Create if necessary.
fid1=fopen([startpath '\codetext.txt'],'at');

marker='**********************************************************************';
line = 0;

% Loop through all mfiles on search path
for i=1:size(files,1)
	filename=files{i};
	fprintf('Compiling %s text.\n',filename);
	% Open mfile(i) and read contents line by line; write lines to codetext.txt
	fid2=fopen(filename,'rt');
	fprintf(fid1,'\n\n%s\n%s\n%s\n\n',marker,filename,marker);
	while line ~=-1
		line=fgets(fid2);
		if line~=-1
			fprintf(fid1,'%s',line);
		end
	end
	line=0;
	% Close mfile(i) and continue
	fclose(fid2);
end
% Close codetext.txt 
fclose(fid1);
fprintf('\n\nDone. Stored in %s.\n\n',[startpath '\codetext.txt']);



function files = filesearch(extension,startpath,recurse)
%
%  Search a specified path (recursively or non-recursively) for all instances of files with a specified extension.
%  Results are returned in a cell array ('files').
%
%  FILES = FILESEARCH(EXTENSION,STARTPATH,1) searches for all instances of *.extension in the starting directory OR
%          in any subdirectory beneath the starting directory.
%  FILES = FILESEARCH(EXTENSION,STARTPATH,0) searches for all instances of *.extension in the starting directory ONLY,
%          and excludes files in subdirectories.
%  Examples:
% 
%  files = filesearch('txt','C:\WinNT',1)
%  mfiles = filesearch('m','D:\Mfiles',0)

% Copyright Brett Shoelson, Ph.D., 5/29/2001 
% brett.shoelson@joslin.harvard.edu

if recurse
	files = {};
	%Create string of recursive directories/subdirectories
	paths = genpath(startpath);
	%Find instances of the path separator
	seplocs = findstr(paths,pathsep);
	%Parse paths into individual (vertically catenated) directories
	if ~isempty(seplocs)
		directories = paths(1:seplocs(1)-1);
		for i = 1:length(seplocs)-1
			directories = strvcat(directories,paths(seplocs(i)+1:seplocs(i+1)-1));
		end
		%Search each directory for instances of *.extension, appending to current list
		for i = 1:size(directories,1)
			%Compile located files as vertically catenated strings
			tmp = dir([deblank(directories(i,:)) '\*.' extension]);
			if ~isempty(tmp),tmp = char(tmp.name);end
			%Update files to reflect newly detected files. (Omit trailing blanks.)
			for j = 1:size(tmp,1)
				files{size(files,1)+1,1} = [deblank(directories(i,:)) '\' deblank(tmp(j,:))];
			end
		end
	end
else %Search non-recursively
	files = {};
	tmp = dir([startpath '\*.' extension]);
	if ~isempty(tmp),tmp = char(tmp.name);end
	for j = 1:size(tmp,1)
		files{size(files,1)+1,1} = [startpath '\' deblank(tmp(j,:))];
	end
end