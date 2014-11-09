function [gaborArray] = makeGaborBank(R, C, u, v)
    Kmax = pi / 2; % max frequency
    f = sqrt( 2 ); % wavelength
    Delt = 2 * pi; % width of gaussian
    Delt2 = Delt * Delt;
    gaborArray = cell(u,v);
    for ix = 1:u
        for j = 1:v
            GW = real(GaborWavelet( R, C, Kmax, f, j, ix, Delt2 )); % Create the Gabor wavelets
            gaborArray{ix,j} = GW;
        end
    end