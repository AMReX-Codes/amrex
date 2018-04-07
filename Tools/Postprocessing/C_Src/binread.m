function [ dim, ngrid, loc, siz, dat] = binread(filename)
%       Read binary AMR data from a file
%
%     dim = number of space dimensions (2 or 3)
%   ngrid = array of number of grids on each level (length=nlev)
%
%  Other quantities are cell arrays of cell arrays:
%   loc{level}{jgrid} = [ xmin xmax ymin ymax ]  ( or [ ... zmin zmax ] )
%   siz{level}{jgrid} = [ nx ny ]                ( or [ nx ny nz ] )
%   dat{ivar}{level}{jgrid} = an array of size nx x ny, or nx x ny x nz
%

if nargin<1
   error('Usage:  binread infile');
end

[ fid, message ] = fopen( filename, 'r' );
if fid == -1
   error(message)
end

%  Read the number of spatial dimensions
[ dim   , count ] = fread(fid,1,'int32');

disp( [ 'dim =' sprintf(' %d',dim) ] );

%  Read the number of levels
[ nlev   , count ] = fread(fid,1,'int32');

%  Read the array sizes
[ ngrid, count ] = fread( fid, nlev, 'int32' );
if count~=nlev;  error('Couldn''t read grid sizes');  end
if any(ngrid<=0)
   error( [ 'Bad number of grids:  ngrid = ' num2str(ngrid) ] );
end

disp( [ 'ngrid =' sprintf(' %d',ngrid) ] );

%  Read the grid physical locations
loc = cell(nlev,1);
for l=1:nlev
   loc{l} = cell(ngrid(l),1);
   for j=1:ngrid(l)
      [ gridloc, count ] = fread( fid, 2*dim, 'double' );
      if count~=2*dim;
         error( [ 'Couldn''t read grid loc at level ' num2str(l) ...
                    '  grid ' num2str(j) ] );
      end
      loc{l}{j} = gridloc;
   end
end

%  Read the grid dimensions
siz = cell(nlev,1);
for l=1:nlev
   siz{l} = cell(ngrid(l),1);
   for j=1:ngrid(l)
      [ gridsiz, count ] = fread( fid, dim, 'int32' );
      if count~=dim;
         error( [ 'Couldn''t read grid dims at level ' num2str(l) ...
                    ' grid ' num2str(j) ] );
      end
      siz{l}{j} = flipud(gridsiz);
      %disp( [ 'l,n,g =' sprintf(' %d %d %d',l,j,siz{l}{j}) ] );
   end
end

%  Read the actual data
nvar = 1;
dat = cell(nvar,1);
for ivar=1:nvar
   dat{ivar} = cell(nlev,1);
   for l=1:nlev
      dat{l} = cell(ngrid(l),1);
      for j=1:ngrid(l)
         ndat = prod(siz{l}{j});
         [ griddat, count ] = fread(fid,ndat,'double');
         if count ~= ndat
            error( [ 'Couldn''t read grid data at level ' num2str(l) ...
                      ' grid ' num2str(j) ] );
         end
         dat{l}{j} = reshape(griddat,siz{l}{j}');
      end
   end
end
