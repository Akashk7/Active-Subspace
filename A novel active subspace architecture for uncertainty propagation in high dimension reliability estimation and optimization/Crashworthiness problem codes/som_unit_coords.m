function Coords = som_unit_coords(topol,lattice,shape)


error(nargchk(1, 3, nargin));

% default values
sTopol = som_set('som_topol','lattice','rect');

% topol
if isstruct(topol), 
  switch topol.type, 
  case 'som_map', sTopol = topol.topol;
  case 'som_topol', sTopol = topol;
  end
elseif iscell(topol), 
  for i=1:length(topol), 
    if isnumeric(topol{i}), sTopol.msize = topol{i}; 
    elseif ischar(topol{i}),  
      switch topol{i}, 
      case {'rect','hexa'}, sTopol.lattice = topol{i}; 
      case {'sheet','cyl','toroid'}, sTopol.shape = topol{i}; 
      end
    end
  end
else
  sTopol.msize = topol;
end
if prod(sTopol.msize)==0, error('Map size is 0.'); end

% lattice
if nargin>1 & ~isempty(lattice) & ~isnan(lattice), sTopol.lattice = lattice; end

% shape 
if nargin>2 & ~isempty(shape) & ~isnan(shape), sTopol.shape = shape; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Action

msize = sTopol.msize;
lattice = sTopol.lattice;
shape = sTopol.shape;

% init variables

if length(msize)==1, msize = [msize 1]; end
munits = prod(msize);
mdim = length(msize);
Coords = zeros(munits,mdim);

% initial coordinates for each map unit ('rect' lattice, 'sheet' shape)
k = [1 cumprod(msize(1:end-1))]; 
inds = [0:(munits-1)]';
for i = mdim:-1:1, 
  Coords(:,i) = floor(inds/k(i)); % these are subscripts in matrix-notation
  inds = rem(inds,k(i)); 
end
% change subscripts to coordinates (move from (ij)-notation to (xy)-notation)
Coords(:,[1 2]) = fliplr(Coords(:,[1 2])); 

% 'hexa' lattice
if strcmp(lattice,'hexa'), 
  % check
  if mdim > 2, 
    error('You can only use hexa lattice with 1- or 2-dimensional maps.');
  end
  % offset x-coordinates of every other row 
  inds_for_row = (cumsum(ones(msize(2),1))-1)*msize(1); 
  for i=2:2:msize(1), 
    Coords(i+inds_for_row,1) = Coords(i+inds_for_row,1) + 1; 
  end
end

% shapes
switch shape, 
case 'sheet', 
  if strcmp(lattice,'hexa'), 
    % this correction is made to make distances to all 
    % neighboring units equal
    Coords(:,2) = Coords(:,2)*sqrt(0.75); 
  end

case 'cyl', 
  % to make cylinder the coordinates must lie in 3D space, at least
  if mdim<3, Coords = [Coords ones(munits,1)]; mdim = 3; end

  % Bend the coordinates to a circle in the plane formed by x- and 
  % and z-axis. Notice that the angle to which the last coordinates
  % are bended is _not_ 360 degrees, because that would be equal to 
  % the angle of the first coordinates (0 degrees).

  Coords(:,1)     = Coords(:,1)/max(Coords(:,1));
  Coords(:,1)     = 2*pi * Coords(:,1) * msize(2)/(msize(2)+1);
  Coords(:,[1 3]) = [cos(Coords(:,1)) sin(Coords(:,1))];
                    
case 'toroid', 

  % NOTE: if lattice is 'hexa', the msize(1) should be even, otherwise 
  % the bending the upper and lower edges of the map do not match 
  % to each other
  if strcmp(lattice,'hexa') & rem(msize(1),2)==1, 
    warning('Map size along y-coordinate is not even.');
  end

  % to make toroid the coordinates must lie in 3D space, at least
  if mdim<3, Coords = [Coords ones(munits,1)]; mdim = 3; end

  % First bend the coordinates to a circle in the plane formed
  % by x- and z-axis. Then bend in the plane formed by y- and
  % z-axis. (See also the notes in 'cyl').

  Coords(:,1)     = Coords(:,1)/max(Coords(:,1));
  Coords(:,1)     = 2*pi * Coords(:,1) * msize(2)/(msize(2)+1);
  Coords(:,[1 3]) = [cos(Coords(:,1)) sin(Coords(:,1))];

  Coords(:,2)     = Coords(:,2)/max(Coords(:,2));
  Coords(:,2)     = 2*pi * Coords(:,2) * msize(1)/(msize(1)+1);
  Coords(:,3)     = Coords(:,3) - min(Coords(:,3)) + 1;
  Coords(:,[2 3]) = Coords(:,[3 3]) .* [cos(Coords(:,2)) sin(Coords(:,2))];
  
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% subfunctions

function C = bend(cx,cy,angle,xishexa)

  dx = max(cx) - min(cx);
  if dx ~= 0, 
    % in case of hexagonal lattice it must be taken into account that
    % coordinates of every second row are +0.5 off to the right
    if xishexa, dx = dx-0.5; end
    cx = angle*(cx - min(cx))/dx; 
  end    
  C(:,1) = (cy - min(cy)+1) .* cos(cx);
  C(:,2) = (cy - min(cy)+1) .* sin(cx);

% end of bend

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

