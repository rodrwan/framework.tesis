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

% u,v,R
gaborArray = makeGaborBank2(7, 3, 9);
% k = 1;
% figure;
% for j = 1:4
%   for i = 1:8
%     gabor_r = gaborArray{i,j};
%     subplot(4, 8, k), subimage(mat2gray(real(gabor_r)));
%     k = k + 1;
%   end
% end
% feat = gabor_conv(scaled, gaborArray);
feat = newFeatures(scaled, 8, 'gabor', 'no');
% imshow(mat2gray(feat(:,:,1)));
k = 1;
figure;
for j = 1:4
  for i = 1:8
    subplot(4,8,k), subimage(mat2gray(feat(:,:,k)));
    k = k + 1;
  end
end