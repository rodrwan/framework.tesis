% load('/Users/rodrwan/Documents/Algorithms/tesis/voc-release5/framework/2011/bottle/bottle_final');
% b_path = '/Users/rodrwan/Documents/Algorithms/tesis/voc-release5/framework/VOC2011/VOCdevkit/VOC2011/JPEGImages/';
% img_id = '2008_000034.jpg';
% img = imread(strcat(b_path, img_id));
%
% im = img);
%
% bbox = process(im, model, -0.5);
%
% showboxes(im, bbox);
%
function test(cls)
% cls = 'bottle';
model_path = '/Users/rodrwan/Documents/Algorithms/tesis/voc-release5/framework/2011/';
model_path = strcat(model_path, '/', cls, '_final.mat')

file = strcat(cls, '_val.txt');
tf_path = '/Users/rodrwan/Documents/Algorithms/tesis/voc-release5/framework/VOC2011/VOCdevkit/VOC2011/ImageSets/Main/';
file_path = strcat(tf_path, file);
[idx, cls] = textread(file_path, '%s %f', 5823);
idc = find(cls == 1);

conf = voc_config();
algorithm = conf.algorithm;
verbose = conf.verbose;
load(model_path);

thresh = -0.5;
for i=1:20;
  b_path = '/Users/rodrwan/Documents/Algorithms/tesis/voc-release5/framework/VOC2011/VOCdevkit/VOC2011/JPEGImages/';
  img_id = strcat(idx{idc(i)}, '.jpg');
  fprintf ('%d.- Testing image: %s\n', i, idx{idc(i)});
  img = imread(strcat(b_path, img_id));
  [ds, bs] = imgdetect(img, model, thresh, algorithm, verbose);
  bbox = process(img, model, thresh, algorithm, verbose);
  showboxes(img, bbox);
  disp('done!');
  pause;
end
