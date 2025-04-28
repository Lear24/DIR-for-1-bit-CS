function beta=assign_initialization(m,d)
% [beta]=initialization(k,d)
% Initialization is performed in a manner that generates 0-1 at each entry.


beta=binornd(1,0.5,m,d);
% beta=binornd(1,0.2,m,d);


end