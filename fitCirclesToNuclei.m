function cMask = fitCirclesToNuclei(bw)

% figure; imagesc(bw);
bou = bwboundaries(bw, 8, 'noholes');
stats = regionprops(~~bw, 'EquivDiameter', 'SubarrayIdx', 'Image');
averageEquivRadius = median([stats.EquivDiameter]/2);
borderThresh = averageEquivRadius; %if an object is within borderThresh px of the image edge, let's not fit a circle to it and leave it be

border = borderImage(bw);
borderDist = bwdist(border);

edgeMask = zeros(size(bw, 1), size(bw, 2));
% lab = bwlabel(bw);

ellipseFrame = zeros(numel(bou), 3);

for i = 1:numel(bou)
%     hold on
    xs = bou{i}(:, 1);
    ys = bou{i}(:, 2);
%     plot(xs, ys, 'g.');
    [xfit,yfit, Rfit]= circfit(xs,ys); %x center y center R 
    
    if borderDist(round(abs(yfit)), round(abs(xfit))) > borderThresh

        ellipseFrame(i, 2) = xfit;
        ellipseFrame(i, 1) = yfit;
        ellipseFrame(i, 3) = Rfit; 
        
    else
        
        reg = stats(i).Image;
%         figure(87); imagesc(reg);
        reg = reg(:);
%         sub = [stats(i).PixelList];
%          subb = stats{i}SubarrayIdx
%         for k = 1:size(reg)
%             if reg(k) == 1
%                 edgeMask(sub(k, 2), sub(k, 1)) = 1;
%             end
% %             figure(88); imagesc(edgeMask);
%         end
        a = [stats(i).SubarrayIdx]; 
        sz = size(stats(i).Image);
        for xx = 1:numel(a{2})
            for yy = 1:numel(a{1})
                edgeMask(a{1}(yy), a{2}(xx)) = reg(sub2ind(sz,yy, xx));
            end
        end
        
    end
    
      
    
    %make circle with center
%     rectangle('position',[xfit-Rfit,yfit-Rfit,Rfit*2,Rfit*2],...
%     'curvature',[1,1],'linestyle','-','edgecolor','r');
end

cMask = makeNuclearMask(ellipseFrame, [size(bw, 1), size(bw, 2)], 'radScale', 1);

cMask = cMask + edgeMask;
% figure; imshowpair(bw, cMask, 'montage');