function [ds, bs, trees] = imgdetect(im, model, thresh, algorithm, verbose)
% Wrapper around gdetect.m that computes detections in an image.
%   [ds, bs, trees] = imgdetect(im, model, thresh)
%
% Return values (see gdetect.m)
%
% Arguments
%   im        Input image
%   model     Model to use for detection
%   thresh    Detection threshold (scores must be > thresh)

% conf = voc_config();
% algorithm = conf.algorithm;
im = color(im);
pyra = featpyramid(im, model, algorithm, verbose);
[ds, bs] = gdetect(pyra, model, thresh);
