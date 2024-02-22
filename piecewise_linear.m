function px = piecewise_linear(x1, x2, x)
	px = x(:);
	for i = 1: size(px, 1)
		if px(i) < x1
			px(i) = 0.;
		elseif(px(i)> x2)
			px(i) = 1.;
		else
			px(i) = (px(i) - x1) / (x2 - x1);
		end
	end
end

