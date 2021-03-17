function [data] = lbl_read2(filename);

%
% function [data] = lbl_read2(filename);
%
% same as lbl_read.m but for single panel files
%
% DCT 11/19/97
%

warning off

h = uicontrol;
set(gcf,'Position',[400 400 350 60],'Color',[0 0 0],...
	'NumberTitle','off','MenuBar','none')
set(h,'Units','normalized','Position',[.05 .1 .9 .8],...
	'BackGroundColor',[.8 .8 .5],'FontSize',16,'FontWeight','b')
set(h,'String','lbl_read.m');drawnow;pause(0.5)

% specify filename if not defined
if nargin == 0;
	[filename,pathname] = uigetfile('*');
	filename = [pathname filename];
end

% open the file as read only and binary format
fid = fopen(filename,'rb');
data.filename = filename;clear filename
if nargin == 2
set(h,'String',['file: ' '' data.filename '' ' opened']);drawnow;pause(0.5)
end

%-----------------------------------------------------------------------
%  read in file_header
%-----------------------------------------------------------------------


set(h,'String','reading file header');drawnow;pause(0.5)

% read in 4 bytes before and after every record 
junk = fread(fid,4,'char');
data.file_header.user_id = setstr(fread(fid,80,'uchar'))';
data.file_header.secant = fread(fid,1,'float64');
data.file_header.p_ave = fread(fid,1,'float32');
data.file_header.t_ave = fread(fid,1,'float32');
molecule_id = zeros(64,8);
for i = 1:64;molecule_id(i,:) = fread(fid,8,'char')';end
data.file_header.molecule_id = setstr(molecule_id);
data.file_header.mol_col_dens = fread(fid,64,'float32');
data.file_header.broad_dens = fread(fid,1,'float32');
data.file_header.dv = fread(fid,1,'float32');
data.file_header.v1 = fread(fid,1,'float64');
data.file_header.v2 = fread(fid,1,'float64');
data.file_header.t_bound = fread(fid,1,'float32');
data.file_header.emis_bound = fread(fid,1,'float32');
data.file_header.LBL_id.hirac = fread(fid,1,'int');
data.file_header.LBL_id.lblf4 = fread(fid,1,'int');
data.file_header.LBL_id.xscnt = fread(fid,1,'int');
data.file_header.LBL_id.aersl = fread(fid,1,'int');
data.file_header.LBL_id.emit = fread(fid,1,'int');
data.file_header.LBL_id.scan = fread(fid,1,'int');
data.file_header.LBL_id.plot = fread(fid,1,'int');
data.file_header.LBL_id.path = fread(fid,1,'int');
data.file_header.LBL_id.jrad = fread(fid,1,'int');
data.file_header.LBL_id.test = fread(fid,1,'int');
data.file_header.LBL_id.merge = fread(fid,1,'int');
data.file_header.LBL_id.scnid = fread(fid,1,'float32');
data.file_header.LBL_id.hwhm = fread(fid,1,'float32');
data.file_header.LBL_id.idabs = fread(fid,1,'int');
data.file_header.LBL_id.atm = fread(fid,1,'int');
data.file_header.LBL_id.layr1 = fread(fid,1,'int');
data.file_header.LBL_id.nlayr = fread(fid,1,'int');
data.file_header.n_mol  = fread(fid,1,'int');
data.file_header.layer = fread(fid,1,'int');
data.file_header.yi1 = fread(fid,1,'float32');
yid = zeros(10,8);
for i = 1:10;yid(i,:) = fread(fid,8,'char')';end
data.file_header.yid = setstr(yid(1:7,:));
% read in 4 bytes before and after every record 
junk = fread(fid,4,'char');
clear molecule_id yid i

% set file_type to "Single"
data.file_type = 'SINGLE';

% Estimate number of spectral points and initialize arrays
n_pts_estimate=ceil((data.file_header.v2-data.file_header.v1)/data.file_header.dv +1.5);
optdepth = zeros(n_pts_estimate,1);

%-----------------------------------------------------------------------
%  read data panel by panel
%-----------------------------------------------------------------------

% initialize counters
data.n_panels = 0;
data.n_pts_read = 0;

% While not end-of-file, read the next panel

end_of_file = 0;
while end_of_file == 0

  % check for end of file
  end_of_file = feof(fid);

  % --------------------------------------------------------
  % read in panel header
  % --------------------------------------------------------

  % read in 4 bytes before and after every record 
  junk = fread(fid,4,'char');

  data.panel_header.v1{data.n_panels+1} = fread(fid,1,'float64');
  data.panel_header.v2{data.n_panels+1} = fread(fid,1,'float64');
  data.panel_header.dv{data.n_panels+1} = fread(fid,1,'float32');
  data.panel_header.n_pts{data.n_panels+1} = fread(fid,1,'int');

  % read in 4 bytes before and after every record 
  junk = fread(fid,4,'char');

  % --------------------------------------------------------
  % read current panel data
  % --------------------------------------------------------

  if data.panel_header.n_pts{data.n_panels+1} ~= -99
	switch upper(data.file_type)
		case 'SINGLE',

		if fix((data.n_panels+1)/50)-(data.n_panels+1)/50 == 0
		set(h,'String',['reading panel # ' num2str(data.n_panels+1) ]);drawnow
		end

		junk = fread(fid,4,'char');
		od = fread(fid,data.panel_header.n_pts{data.n_panels+1},'float32');
		junk = fread(fid,4,'char');

		pt1 = data.n_pts_read + 1;
		pt2 = data.n_pts_read + data.panel_header.n_pts{data.n_panels+1};
		optdepth(pt1:pt2) = od;
		clear od

		% increment number of points read
		data.n_pts_read = data.n_pts_read + ...
			data.panel_header.n_pts{data.n_panels+1};

		% increment panel number
		data.n_panels = data.n_panels +1;

	end
  end
end  % on end-of-file while loop

set(h,'String','close file and return data');drawnow

data.od = optdepth(1:data.n_pts_read);
data.v1 = data.panel_header.v1{1};
data.v2 = data.panel_header.v2{data.n_panels};
data.v = linspace(data.v1,data.v2,data.n_pts_read)';

fclose(fid);

close(gcf)
drawnow
