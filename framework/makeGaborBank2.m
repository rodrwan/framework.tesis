function GC = makeGaborBank2(u, v, R)
  Kmax = pi/2;  % max frequency
  f = sqrt(2);  % wavelength
  U = 0:u;      % orientation
  V = 0:v;      % scales
  % this sigma is dynamic couse we need to changes the filter size.
  % R = 8*n*pi;
  % => n = R/(8*pi);
  % => sigma = n*pi
  sigma = R/(8); % width of gaussian
  %% Defaults
  loadFilter = true;
  filterFile = 'gabor.mat';
% global loadFilter; global filterFile;
% Kmax: max frequency (default: pi/2)
% f: is the wavelet value (default: sqrt(2))
% u: orientation
% v: frequency
% sigma: the width of the gaussian

% if the set exist load it
if loadFilter && exist(filterFile,'file')
  load(filterFile)
end
if ~exist('GC','var')
  R=fix(8*sigma);
  Delt2 = sigma * sigma;

  n_u = length(U);
  n_v = length(V);
  GC = cell(n_u, n_v);
  mesh = -R/2 + 1 : R/2;
  for i = 1:n_u
    for j = 1:n_v
      GW = zeros (R,R);
      u = U(i);
      v = V(j);
      k = ( Kmax/(f^v) )*exp( 1i*u*pi/8 );% Wave Vector
      kn2 = ( abs(k) )^2;
      for m = mesh
          for n = mesh
              GW(m+R/2,n+R/2) = (kn2/Delt2)*exp(-0.5*kn2*(m^2 + n^2)/Delt2)*(exp(1i*(real(k)*m + imag(k)*n))-exp(-0.5*Delt2));
          end
      end
      GC{i, j} = GW;
    end
  end
  % save the filter to a file to speed up
  if loadFilter
    save(filterFile, 'GC');
  end
end
end