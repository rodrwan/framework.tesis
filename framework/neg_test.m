clear all;
close all;
% compile_f
% function [feat] = neg_test(algorithm)
img = '/Users/rodrwan/Documents/Algorithms/tesis/voc-release5/framework/VOC2011/VOCdevkit/VOC2011/JPEGImages/2008_000008.jpg';
% img = '/Users/rodrwan/Documents/Algorithms/tesis/voc-release5/framework/VOC2011/VOCdevkit/VOC2011/JPEGImages/2008_000162.jpg';

imr = imread(img);
imshow(imr);
fprintf('Image dimention: %dx%d\n', size(imr, 1), size(imr, 2));
scaled = double(imr);
sbin = 8;

% fprintf('\nPrint gabor result:\n');
verbose = 'yes';
algorithm = 'gabor';
feat = newFeatures(scaled, sbin, algorithm, verbose);
% imshow(mat2gray(feat(:,:,1)));
figure;
for k = 1:32
  subplot(4,8,k), subimage(mat2gray(feat(:,:,k)));
end