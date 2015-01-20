function graph

opt = csvread('optimization.csv', 1, 0);
tar = csvread('Tarakanov.csv', 1, 0);

oEta = opt(:, 1);
oOme = opt(:, 2);
oLam0 = opt(:, 3);
oLam1 = opt(:, 4);
oT = opt(:, 5);

tT = tar(:, 1);
tOme = tar(:, 4);
tLam0 = tar(:, 6);
tLam1 = tar(:, 7);

h1 = figure();
plot(oT, oOme, tT, tOme);

% h2 = figure();
% plot(oT, oOme, oT, oLam0./(oLam0 + oLam1), oT, oLam1./(oLam0 + oLam1));
% 
% h3 = figure();
% plot(oT, oLam0./(oLam0 + oLam1), oT, oLam1./(oLam0 + oLam1));
% 
% h4 = figure();
% plot(oT, oOme, oT, (oLam1 + oLam0), oT, (1 - oOme - oLam0 - oLam1),'.');

h5 = figure();
a1 = oLam1 + oLam0;
a2 = ones(size(oLam1)) - oOme;
a3 = ones(size(oLam1));
plot(oT, a1, oT, a2, oT, a3);

end