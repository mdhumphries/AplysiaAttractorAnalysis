function cmap = sequentialYOB(N)

x = 0:1/N:1;
R = 1 - 0.392*(1 + erf((x - 0.869)/0.255));
G = 1.021 - 0.456*(1 + erf((x - 0.527)/0.376));
B = 1 - 0.493* (1 + erf((x -0.272)/0.309));

cmap = [R' G' B'];
