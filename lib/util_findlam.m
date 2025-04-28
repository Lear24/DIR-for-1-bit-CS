function lam = util_findlam(Accuracy_rate,lam_information,lambda)

ind1 = find(Accuracy_rate == min(Accuracy_rate(1:lam_information(end))));

ind2 = find(lam_information >= ind1(1));

lam = lambda(ind2(1));

end