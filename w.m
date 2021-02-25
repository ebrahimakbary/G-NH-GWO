function y = w(x,c11,c22)

y = zeros(length(x),1);
for k = 1:length(x)
	y(k) = sum(c11 .* cos(c22.*x(:,k)));
end