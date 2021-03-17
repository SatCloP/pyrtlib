function files = filesearch(extension,startpath,recurse)
%
%  Search a specified path (recursively or non-recursively) for all instances of files with a specified extension.
%  Full file names are returned in a cell array ('files').
%
%  FILES = FILESEARCH(EXTENSION,STARTPATH,1) searches for all instances of *.extension in the starting directory OR
%          in any subdirectory beneath the starting directory.
%  FILES = FILESEARCH(EXTENSION,STARTPATH,0) searches for all instances of *.extension in the starting directory ONLY,
%          and excludes files in subdirectories.
%  Examples:
% 
%  files = filesearch('txt','C:\WinNT',1)
%  mfiles = filesearch('m','D:\Mfiles',0)
%
% Copyright Brett Shoelson, Ph.D. 
% Shoelson Consulting
% brett.shoelson@joslin.harvard.edu
% 5/29/2001

if nargin ~= 3
	error('Requires 3 input arguments');
end

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
else
	files = {};
	tmp = dir([startpath '\*.' extension]);
	if ~isempty(tmp),tmp = char(tmp.name);end
	for j = 1:size(tmp,1)
		files{size(files,1)+1,1} = [startpath '\' deblank(tmp(j,:))];
	end
end