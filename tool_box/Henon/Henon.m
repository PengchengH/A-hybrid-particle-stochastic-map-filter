function OUT = Henon(IN)
OUT = zeros(size(IN));
OUT(1,:) = 1-1.4*IN(1,:).^2+IN(2,:);
OUT(2,:) = 0.3*IN(1,:);