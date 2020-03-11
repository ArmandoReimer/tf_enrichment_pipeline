sigma = 7;
kmask = imsegkmeans(single(imgaussfilt(im,sigma)),3)==3;
snakesFun = @(b, s, sigma) activecontour(imgaussfilt(im, sigma), kmask, 'Chan-Vese', 'ContractionBias', b, 'SmoothFactor', s);
b = -.4; s = .1;
kMaskRefined = snakesFun(b, s, sigma/2);
refMask = wshed(kMaskRefined);



%just to see what this looks like
% immask = ~kmask.*im

%use the mask as a starting point and grow it outward to better capture the
%nuclear boundaries

%potentially i could create a loop over smooth factors and check if they
%change the number of regions in the mask. then i pick the highest
%smoothing that retains the number of regions or something similar.
% figure; imagesc(kMaskRefined); colorbar;
% figure; imshowpair(kmask, kMaskRefined)

% immaskRef = ~kMaskRefined.*im;
% figure; imagesc(immaskRef); colorbar;
% figure; imshowpair(immaskRef, im, 'Scaling', 'joint')
% figure; imshowpair(immaskRef, im, 'montage','Scaling', 'joint')
b = wshed(kMaskRefined);





% g = imgradient(imgaussfilt(immask,8));
% t = graythresh(g);
% bw = imbinarize(g);
% figure; imagesc(bw);
% %nah. try drawing a larger bounding box and finding the centroid (maybe
% %smooth before centroiding)