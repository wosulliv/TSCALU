x = 4096;
y = 4096;
n = 4096;

C = gallery('circul',x);
imagesc(C);
title('Circulant Matrix');
colorbar
cC = cond(C);

D = gallery('dorr',x);
imagesc(D);
title('Diagonally Dominant Matrix');
colorbar
cD = cond(D);

J = gallery('jordbloc',x);
imagesc(J);
title('Jordan Block Matrix')
colorbar
cJ = cond(J);

K = gallery('kms',x);
imagesc(K);
title('Kac-Murdock-Szegö Matrix')
colorbar
cK = cond(K);

N = gallery('neumann',n);
imagesc(N);
title('Neumann Matrix');
colorbar
cN = cond(N);

RGM = normrnd(0,1,16384);
imagesc(RGM);
title('Random Gaussian Matrix');
colorbar
cRGM = cond(RGM);