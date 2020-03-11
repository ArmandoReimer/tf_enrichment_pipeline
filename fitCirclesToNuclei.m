function cMask = fitCirclesToNuclei(mask)

% figure; imagesc(bw);
xDim = size(mask, 1);
yDim = size(mask, 1);

boundaryCell = bwboundaries(mask, 8, 'noholes');
stats = regionprops(~~mask, 'EquivDiameter', 'SubarrayIdx', 'Image');
averageEquivRadius = median([stats.EquivDiameter]/2);
%if an object is within borderThresh px of the image edge, let's not fit a circle to it and leave it be
borderThresh = averageEquivRadius; 
border = borderImage(mask);
borderDist = bwdist(border);

edgeMask = false(xDim, yDim);

ellipseFrame = zeros(numel(boundaryCell), 3);

for i = 1:numel(boundaryCell)
    %     hold on
    xs = boundaryCell{i}(:, 1);
    ys = boundaryCell{i}(:, 2);
    %     plot(xs, ys, 'g.');
    
    %[x center y center R]
    [xfit,yfit, Rfit]= circfit(xs,ys);
    
    xSub = max(round(abs(xfit)), xDim);
    ySub = max(round(abs(yfit)), yDim);
    
    isNearBorder =  borderDist(ySub, xSub) > borderThresh;
    
    if isNearBorder
        
        ellipseFrame(i, 2) = xfit;
        ellipseFrame(i, 1) = yfit;
        ellipseFrame(i, 3) = Rfit;
        
    else
        
        reg = stats(i).Image;
        %         figure(87); imagesc(reg);
        reg = reg(:);
        % %             figure(88); imagesc(edgeMask);
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

cMask = makeNuclearMask(ellipseFrame, [size(mask, 1), size(mask, 2)], 'radiusScale', 1);

cMask = cMask + edgeMask;
% figure; imshowpair(bw, cMask, 'montage');