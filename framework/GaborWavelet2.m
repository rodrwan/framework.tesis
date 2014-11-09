function GC = makeGaborBank2(Kmax, f, U, V, sigma)
  %% Defaults
  global loadFilter; loadFilter = true;
  global filterFile; filterFile = 'gabor.mat';
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
if ~exist('GW','var')
  R=fix(8*sigma);
  Delt2 = sigma * sigma;

  n_u = length(U);
  n_v = length(V);

  GW = zeros (R,R,n_u,n_v);
  mesh = -R/2 + 1 : R/2;

  for i = 1:n_u
    for j = 1:n_v
      u = U(i);
      v = V(j);
      k = ( Kmax/(f^v) )*exp( 1i*u*pi/8 );% Wave Vector
      kn2 = ( abs(k) )^2;
      for m = mesh
          for n = mesh
              GW(m+R/2,n+R/2) = (kn2/Delt2)*exp(-0.5*kn2*(m^2 + n^2)/Delt2)*(exp(1i*(real(k)*m + imag(k)*n))-exp(-0.5*Delt2));
          end
      end
      GC{i, j} = GW
    end
  end
  % save the filter to a file to speed up
  if loadFilter
    save(filterFile,'GW');
  end
end
end