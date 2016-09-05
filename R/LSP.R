#' Lomb Scargle Periodogram
#'
#' Calculates the Lomb Scargle Periodogram using the algorithm described by K. Hocke and N. Kämpfer (2009).
#' @param t Numeric vector of timepoints.
#' @param y Numeric vector of values corresponding to t.
#' @param ofac Oversampling factor. Recommend 2-4. Higher improves granularity, but may cause artefacts.
#' @param MCsim MCsim object returned by function "monte_carlo_lsp". If not NULL will return p-values.
#' @param zero_factor Integer, zero padding factor. Pads fourier series to increase resolution of approximated function. Number of zeros is equal to zero_factor*length of (unpadded) fourier series. Value of 0 means no padding.
#' @param pval_thresh Filter fourier series based on peak strength p value. This will affect the reconstructed series.
#' @return A list with the following elements:\cr
#' t - original t.\cr
#' y - original y.\cr
#' t2 - approximated t.\cr
#' y2 - reconstructed y. plot against t2.\cr
#' frequency - frequency range used.\cr
#' spectral_power_density - Plot against Frequency for the periodogram\cr
#' peak_info - data.frame with index, frequency, spd, and (optional) pvalue of each detected peak.\cr
#' peak_spd - highest spd\cr
#' peak_period - period corresponding to peak_spd\cr
#' fourier - fourier series\cr
#' fourier_freq frequencies for fourier series\cr
#' @references Hocke, K., and N. Kämpfer. "Gap filling and noise reduction of unevenly sampled data by means of the Lomb-Scargle periodogram." Atmospheric Chemistry and Physics 9.12 (2009): 4197-4206.
#' @export
lsp = function(t,y,ofac=2,MCsim=NULL,zero_factor=0,pval_thresh=1) {
  xstart=t[1];
  x=t-xstart;

  hifac=1;
  twopi=2*pi;
  n=length(x);
  nout=0.5*ofac*hifac*n;
  nmax=nout;

  wi=rep(0,nmax);
  wpi=wi;
  wpr=wi;
  wr=wi;
  wtemp=wi;
  px=wi;
  py=wi;
  ph=wi;
  ph1=wi;
  Fx=wi;
  Fy=wi;
  ave=mean(y);
  vari=var(y);

  xmax=max(x);
  xmin=min(x);
  xdif=xmax-xmin;
  xave=0.5*(xmax+xmin);

  pymax=0.;
  pnow=1/(xdif*ofac);
  arg=twopi*((x-xave)*pnow);
  wpr=-2*sin(0.5*arg)^2;
  wpi=sin(arg);
  wr=cos(arg);
  wi=wpi;
  yy=(y-ave);
  for(i in 1:nout) {
    px[i]=pnow;
    sumsh=sum(wr*wi);
    sumc=sum((wr-wi)*(wr+wi));
    wtau=0.5*atan2(2.*sumsh,sumc);
    swtau=sin(wtau);
    cwtau=cos(wtau);
    ss=wi*cwtau-wr*swtau;
    cc=wr*cwtau+wi*swtau;
    sums=sum(ss^2);
    sumc=sum(cc^2);
    sumsy=sum(yy*ss);
    sumcy=sum(yy*cc);
    wtemp=wr;
    wr=wr*wpr-wi*wpi+wr;
    wi=wi*wpr+wtemp*wpi+wi;
    iy=sumsy/sqrt(sums); #% imaginary part of Lomb-Scargle spectral component
    ry=sumcy/sqrt(sumc); #% real part
    py[i]=0.5*(ry^2+iy^2)/vari; #% power
    #% here, the FFT phase is computed from the Lomb-Scargle Phase
    #% at each new frequency 'pnow' by adding the phase shift 'arg0'
    phLS=atan2(iy,ry);            #% phase of Lomb-Scargle spectrum
    arg0=twopi*(xave+xstart)*pnow +wtau;#  % phase shift with respect to 0
    arg1=twopi*xave*pnow +wtau;   #% phase shift for FFT reconstruction
    ph[i]=(phLS+arg0) %% twopi;  #% phase with respect to 0
    ph1[i]=(phLS+arg1) %% twopi; #% phase for complex FFT spectrum
    pnow=pnow+1./(ofac*xdif);    # next frequency
  }

  dim=2*nout+1;    #%dimension of FFT spectrum
  fac=sqrt(vari*dim/2);
  a=fac*sqrt(py);    #% amplitude vector for FFT
  Fx=a*cos(ph1); #% real part of FFT spectrum
  Fy=a*sin(ph1); #% imaginary part of FFT spectrum
  ph=ph + 5*twopi %% twopi;    #% for value range 0,..., 2 pi
  wk1=px ;
  wk2=py;


  Fxr=rev(Fx)[-1]
  Fyr=rev(Fy)[-1]
  #%complex Fourier spectrum which corresponds to the Lomb-Scargle periodogram:
  fou = c(complex(real=ave,imaginary=0),complex(real=Fx,imaginary=-Fy),complex(real=Fxr,imaginary=Fyr));

  FFreq = c(0,wk1[-length(wk1)]);

  PeakIndex = match(max(wk2), wk2)
  PeakPeriod <- 1 / wk1[PeakIndex]

  N = length(t);
  NF = length(fou)
  pp = pracma::findpeaks(wk2,threshold = 0)


  if(is.null(MCsim)) {
    peak_info = data.frame(idx=pp[,2],freq=wk1[pp[,2]],spd=pp[,1]);
  } else {
    spd = pp[,1];
    pvals = numeric(length(spd));
    for(i in 1:length(spd)) {
      pvals[i] = calc_p(MCsim,spd[i]);
    }
    peak_info = data.frame(idx=pp[,2],freq=wk1[pp[,2]],spd=pp[,1],pval=pvals);
  }
  peak_info = peak_info[order(peak_info[,3],decreasing = T),]



  idxs = peak_info[["idx"]]+1;

  if(!is.null(MCsim)) {
    idxs = idxs[peak_info[["pval"]] <= pval_thresh]
    fou[-idxs]  = complex(1,0,0);
    fou[(NF/2+1):NF] = rev(Conj(fou[2:(NF/2+1)]));
  }

  if(length(idxs)==0){
    significant_peaks_found = F;
  } else {
    significant_peaks_found = T;
  }


  #zero-pad
  fourier = c(fou[1:(NF/2)],rep(complex(1,0,0),NF*zero_factor),fou[(NF/2+1):NF])

  if(significant_peaks_found) {
    #ifft
    inverse_F = pracma::ifft(fourier);

    y2 = Re(inverse_F)[1:(length(inverse_F)/ofac)]
    y2 = y2 + Re(fourier[1])
    # y2[length(y2)] = y2[1];
    t2 = seq(t[1],t[length(t)], length.out = length(y2)+1);
    t2 = t2[-length(t2)];
    y2i = approx(t2,y=y2,xout=t);

    rrr = MASS::rlm(y ~ y2i$y)
    y2 = rrr$coefficients[2]*y2 + rrr$coefficients[1];
  } else {
    t2 = numeric(0)
    y2 = numeric(0)
  }


  return( list( t=t,
                y=y,
                t2=t2,
                y2=y2,
                frequency=wk1,
                spectral_power_density=wk2,
                peak_info=peak_info,
                peak_spd=wk2[PeakIndex],
                peak_period=PeakPeriod,
                peak_freq=wk1[PeakIndex],
                fourier = fourier,
                fourier_freq = FFreq

  ) )
}

lsp_mc = function(t,N_runs,ofac) {
  mxs = rep(0,N_runs);
  # pb = txtProgressBar(min = 0, max = N_runs, initial = 0, char = "=", width = 50, style = 3)
  for(j in 1:N_runs) {

    y=rnorm(length(t),mean=0,sd=1)
    xstart=t[1];
    x=t-xstart;

    hifac=1;
    twopi=2*pi;
    n=length(x);
    nout=0.5*ofac*hifac*n;
    nmax=nout;

    wi=rep(0,nmax);
    wpi=wi;
    wpr=wi;
    wr=wi;
    wtemp=wi;
    px=wi;
    py=wi;
    ph=wi;
    ph1=wi;
    Fx=wi;
    Fy=wi;
    ave=mean(y);
    vari=var(y);

    xmax=max(x);
    xmin=min(x);
    xdif=xmax-xmin;
    xave=0.5*(xmax+xmin);

    pymax=0.;
    pnow=1/(xdif*ofac);
    arg=twopi*((x-xave)*pnow);
    wpr=-2*sin(0.5*arg)^2;
    wpi=sin(arg);
    wr=cos(arg);
    wi=wpi;
    yy=(y-ave);
    for(i in 1:nout) {
      px[i]=pnow;
      sumsh=sum(wr*wi);
      sumc=sum((wr-wi)*(wr+wi));
      wtau=0.5*atan2(2.*sumsh,sumc);
      swtau=sin(wtau);
      cwtau=cos(wtau);
      ss=wi*cwtau-wr*swtau;
      cc=wr*cwtau+wi*swtau;
      sums=sum(ss^2);
      sumc=sum(cc^2);
      sumsy=sum(yy*ss);
      sumcy=sum(yy*cc);
      wtemp=wr;
      wr=wr*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
      iy=sumsy/sqrt(sums); #% imaginary part of Lomb-Scargle spectral component
      ry=sumcy/sqrt(sumc); #% real part
      py[i]=0.5*(ry^2+iy^2)/vari; #% power

      pnow=pnow+1./(ofac*xdif);    # next frequency
    }

    # print(j)
    mxs[j] = max(py);
    # setTxtProgressBar(pb, value=j)
  }

  return(mxs)
}

#' LSP + Monte Carlo
#'
#' Runs the LSP algorithm a number of iterations using random values on the specified timepoints. Records highest peak at each iteration.
#' @param t Numeric vector of timepoints
#' @param Nruns number of monte carlo iterations
#' @param ofac oversampling factor, should be same as in lsp.
#' @param N_cores number of cores to use. Default 1.
#' @return Returns raw data of the monte carlo runs. To be passed on to the "MCsim" parameter of function "lsp".
#' @export
monte_carlo_lsp = function(t,Nruns,ofac,N_cores = 1) {
  runs_per_core = floor(Nruns/N_cores);
  remainder = Nruns %% N_cores;
  runs = rep(runs_per_core,N_cores);

  if(remainder>0) {
    for(i in 1:remainder) {
      runs[i] = runs[i]+1;
    }
  }

  cl <- makeCluster(N_cores);
  result = unlist(parSapply(cl,runs,lsp_mc,t=t,ofac=ofac))
  stopCluster(cl)

  return(list(t=t,Nruns=Nruns,ofac=ofac,peaks=result));
}

#' Calculate P-value
#'
#' Calculates p_value using monte carlo simulation data and a sample spectral power density value.
#' @param MCsim output from monte_carlo_lsp
#' @param test_spd spd to convert to p-value
#' @export
calc_p = function(MCsim,test_spd) {
  return( sum(MCsim$peaks>test_spd)/MCsim$Nruns )
}

NHorneBaliunas <- function(N)
{
  Nindependent <- trunc(-6.362 + 1.193*N + 0.00098*N^2)
  if (Nindependent < 1)
  {
    Nindependent <- 1   # Kludge for now
  }
  return(Nindependent)
}
