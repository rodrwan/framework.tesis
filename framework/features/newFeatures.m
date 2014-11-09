function new_feature = newFeatures(scaled, sbin, algorithm, verbose)
    if strcmp(algorithm, 'gabor') == 1
        gaborArray = makeGaborBank2(7, 3, 9);
        new_feature = features(scaled, sbin, gaborArray, 'gabor', verbose);
    else
        gaborArray = [];
        new_feature = features(scaled, sbin, gaborArray, 'hog', verbose);
    end