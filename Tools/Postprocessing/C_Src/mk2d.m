
% Erase all earlier variables from matlab.
clear all;

% Clear the plotting window.
clf;

% This is the (integer) time stamp on the .mat file.
i = 50;

% Change "density" to the appropriate variable
%fname = sprintf('plt%04d_density.mat',i);
fname = 'plt03310_temp.mat';

% Read in the number of grids at each level, the locations and sizes 
%  of those grids, and the data on each grid.
[dim ngrid loc siz dat] = binread(fname);

% Number of levels
nlev = size(ngrid);

% The size of the graph will be the size of the first level 0 grid.
% If there is more than one level 0 grid, this will need to be modified.
xlo = loc{1}{1}(1);
xhi = loc{1}{1}(2);
ylo = loc{1}{1}(3);
yhi = loc{1}{1}(4);

for n=1:ngrid(1)
  xlo = min(xlo,loc{1}{n}(1));
  xhi = max(xhi,loc{1}{n}(2));
  ylo = min(ylo,loc{1}{n}(3));
  yhi = max(yhi,loc{1}{n}(4));
end

axis([xlo xhi ylo yhi]);

ncont = 10;
showBoxes = 1;

zmin = min(dat{1}{1}(:));
zmax = max(dat{1}{1}(:));
for l=1:nlev
  num_grids = ngrid(l);
  for n=1:num_grids
    zmin = min(zmin,min(dat{l}{n}(:)));
    zmax = max(zmax,max(dat{l}{n}(:)));
  end
end
zdelta = (zmax - zmin)/ncont;

hold on;

% Contour the data at each level and draw the grid boxes.
for l=1:nlev
  num_grids = ngrid(l);
  for n=1:num_grids
    xlo = loc{l}{n}(1);
    xhi = loc{l}{n}(2);
    ylo = loc{l}{n}(3);
    yhi = loc{l}{n}(4);

    xvec = [xlo xlo xhi xhi xlo];
    yvec = [ylo yhi yhi ylo ylo];
    if (l~=0)
      fill(xvec,yvec,'w','edgecolor','w');
    end
    
    if (showBoxes~=0)
      plot(xvec,yvec,'k-');
    end
    
  end
end

for l=1:nlev
  num_grids = ngrid(l);
  for n=1:num_grids
    xlo = loc{l}{n}(1);
    xhi = loc{l}{n}(2);
    ylo = loc{l}{n}(3);
    yhi = loc{l}{n}(4);
    nx  = siz{l}{n}(2);
    ny  = siz{l}{n}(1);
    dx = (xhi-xlo)/nx;
    dy = (yhi-ylo)/ny;
    [x,y] = meshgrid(xlo+dx/2:dx:xhi-dx/2, ylo+dy/2:dy:yhi-dy/2);
    z = dat{l}{n};

    contour(x,y,z,zmin:zdelta:zmax);

  end
end

axis image;

