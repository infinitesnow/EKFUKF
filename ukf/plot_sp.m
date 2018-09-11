set(0, 'currentfigure', sp_figure);

for ii = 1:size(Sigma_points,2)
    scatter3(Sigma_points(1,ii),Sigma_points(2,ii),Sigma_points(3,ii))
end
SP(k) = getframe;