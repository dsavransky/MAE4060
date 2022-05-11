function [truncpath, path2circ1, path2circ2, circintinds] = ...
    findLambertIntersections(r1v,r2v,depcirc,path)

%find intersections with endpoints and departure circle
[~,ind1] = min(sqrt((path(1,:) - r1v(1)).^2 + (path(2,:) - r1v(2)).^2));
[~,ind2] = min(sqrt((path(1,:) - r2v(1)).^2 + (path(2,:) - r2v(2)).^2));
[xint,yint,pathind,circind] = intersections(path(1,:),path(2,:),...
    depcirc(1,:), depcirc(2,:));

inds = sort([ind1,ind2]);
truncpath = path(:,inds(1):inds(2));
inds = sort([ind1,round(pathind(1))]);
path2circ1 = path(:,inds(1):inds(2));
inds = sort([ind1,round(pathind(2))]);
path2circ2 = path(:,inds(1):inds(2));
circintinds = round(circind);



