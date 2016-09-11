library(LSP)

#set up testing parameters
ofac = 8;
amp=3;
sd = 1.0;
n = 100; #initial series length
n_out = 50; #sereis length after random sampling
period = 45;
freq = 1/period;
offset = 25;
phase = pi/2;
p_thresh = 0.05;
N_mc_runs = 1000;
zero_factor = 2;
N_cores = 1;

#create testing data
x = seq(0,n-1);
x = x + offset;
y = amp*cos(2*pi*freq*x - phase);

#random sampling to simulate missing data/irregular timepoints
idx = seq(1,n);
idx = sample(idx,n_out);
idx = sort(idx)

#add some noise
y = rnorm(n,0,sd) + y;

#final data
y = y[idx];
x = x[idx];

#monte carlo
mc_sim = monte_carlo_lsp(x,N_mc_runs,ofac,N_cores);

#LSP
# result = lsp(x,y,ofac=ofac,mc_sim = mc_sim,zero_factor = zero_factor, thresh = p_thresh,verbose=1);
result = lsp(x,y,ofac=ofac,mc_sim = mc_sim,zero_factor = zero_factor, thresh = 0.05,verbose=1);

#plot input data vs reconstructed signal
plot(result$t,result$y,col="red")
lines(result$t2,result$y2,col="green")


