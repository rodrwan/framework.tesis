function app(img, model, thresh)
startup;

fprintf('compiling the code...');
compile;
fprintf('done.\n\n');

format shortG
load(model);
im = imread(img);
bbox = process(im, model, thresh)
% showboxes(im, bbox);