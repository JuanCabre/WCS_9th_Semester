function grid = initPropaFriis(gridLength, z, objectLoc, Pt, alpha, lambda, Gr, Gt)
	%Fill propogation into matrix
	[i, k] = ndgrid(1:gridLength(1), 1:gridLength(2));
	%Friis equation
	grid(:,:) = 10*log10((Pt*(lambda./((4*pi*sqrt( ( (objectLoc(1) - i).*z ).^2 + ( (objectLoc(2) - k).*z ).^2 )))).^alpha)*Gr*Gt);
	%Insert full transmit power into target position
	grid(objectLoc(1), objectLoc(2)) = 10*log10(Pt);
end