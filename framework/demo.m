function [p, r] = demo()
startup;

fprintf('compiling the code...');
% compile;
fprintf('done.\n\n');

load('/home/dev/rfuenzalida-tesis/tesis/framework/2011/bottle_final');
model.vis = @() visualizemodel(model, 1:2:length(model.rules{model.start}));
cls = model.class;
clf;
% load and display model
model.vis();
disp([cls ' model visualization']);
disp('press any key to continue'); pause;
disp('continuing...');

tf_path = '/home/dev/rfuenzalida-tesis/tesis/framework/VOC2011/VOCdevkit/VOC2011/ImageSets/Main/bottle_trainval.txt';
[idx, clss] = textread(tf_path, '%s %f', 5823);
idc = find(clss == 1);
dtc = zeros(size(idc,1), 1);
% tiempo total
tt = 0;
for i=1:size(idc,1);
  b_path = '/home/dev/rfuenzalida-tesis/tesis/framework/VOC2011/VOCdevkit/VOC2011/JPEGImages/';
  img_id = strcat(idx{idc(i)}, '.jpg'); % idc(i)
  fprintf ('%d.- Testing image: %s,', i, idx{idc(i)}); %idc(i)
  img = strcat(b_path, img_id);

  % read annotation
  a_path = '/home/dev/rfuenzalida-tesis/tesis/framework/VOC2011/VOCdevkit/VOC2011/Annotations/';
  img_an = strcat(idx{idc(i)}, '.xml'); % idc(i)
  filename = strcat(a_path, img_an);
  try
    annotation = PASreadrecord(filename);
  catch
    error('Failed to read XML file %s.',filename);
  end

  bbox = annotation.objects.bbox;
  [cls, et] = test2(img, model, -0.5, bbox);
  fprintf(' is detected? %d,', cls);
  dtc(i) = cls;
  fprintf(' elapsed time: %.4f secs\n', et);
  tt = tt + et;
end
fprintf('Mean time: %4f secs\n', tt/size(idc,1));

% precision
scores = dtc;
targets = zeros(size(dtc,1), 1);
% recall
for i=1:size(idc,1);
    if dtc(i) == 1 && dtc(i) == clss(idc(i))
        targets(i) = 1;
    else
        targets(i) = -1;
    end
end

figure
[Xpr,Ypr,Tpr,AUCpr] = perfcurve(targets, scores, 1, 'xCrit', 'reca', 'yCrit', 'prec');
plot(Xpr,Ypr)
xlabel('Recall'); ylabel('Precision')
title(['Precision-recall curve (AUC: ' num2str(AUCpr) ')'])
% load('INRIA/inriaperson_final');
% model.vis = @() visualizemodel(model, ...
%                   1:2:length(model.rules{model.start}));
% test('000061.jpg', model, -0.3);
%
% load('VOC2007/person_grammar_final');
% model.class = 'person grammar';
% model.vis = @() visualize_person_grammar_model(model, 6);
% test('000061.jpg', model, -0.6);
%
% load('VOC2007/bicycle_final');
% model.vis = @() visualizemodel(model, ...
%                   1:2:length(model.rules{model.start}));
% test('000084.jpg', model, -0.3);

function test(imname, model, thresh)
    cls = model.class;
    fprintf('///// Running demo for %s /////\n\n', cls);

    % load and display image
    clf;
    im = imread(imname);
    imshow(im);
    % axis equal;
    % axis on;
    disp('input image');
    disp('press any key to continue'); pause;
    disp('continuing...');

    % detect objects
    [ds, bs] = imgdetect(im, model, thresh);

    if ~isempty(ds) || ~isempty(bs)
        top = nms(ds, 0.5);
        clf;
        if model.type == model_types.Grammar
          bs = [ds(:,1:4) bs];
        end
        showboxes(im, reduceboxes(model, bs(top,:)));
        disp('detections');
        disp('press any key to continue'); pause;
        disp('continuing...');

        bbox = process(im, model, -0.5);
        showboxes(im, bbox);

    else
       disp('No se detecto nada');
    end
    fprintf('\n');


function [dtc, th] = test2(imname, model, thresh, bbox)
  warning('off','all');
  x1 = bbox(1);
  y1 = bbox(2);
  x2 = bbox(3);
  y2 = bbox(4);
  conf = voc_config();
  algorithm = conf.algorithm;
  verbose = conf.verbose;

  im = imread(imname);
  th = tic();
  bbox = process(im, model, thresh, algorithm, verbose);
  th = toc(th);

  if ~isempty(bbox)
    r_x1 = bbox(1);
    r_y1 = bbox(2);
    r_x2 = bbox(3);
    r_y2 = bbox(4);
    score = bbox(6);
    if x1*0.5 < r_x1 && x2*0.5 < r_x2 && y1*0.5 < r_y1 && y2*0.5 < r_y2
      dtc = score;
      % clf;
      % imshow(im);
      % showboxes(im, bbox);
    else
      dtc = 0;
    end
  else
    dtc = 0;
  end
  % pause(1);
