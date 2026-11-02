/* Headless adapter to the existing STRF rffit numerical engine.
 * Uses rffit.c and PGPLOT; units internally kHz, exported Hz.
 * argv: data catalog output.tle residuals.csv mode[fixed|ma|mamm|raan|maraan|mammraan] fit_cut_mjd [plot_device]
 */
#define main strf_interactive_main
#include "../rffit.c"
#undef main
extern int STRF_FIT_STATUS;
/* TCA: approaching-to-receding radial velocity zero above the horizon.
 * Same physical definition as interactive rffit, refined by bisection. */
static double pass_tca(site_t site) {
  double prev=d.mjdmin,v0,v,az,el;
  velocity(orb,prev,site,&v0,&az,&el);
  for(int j=1;j<=NMAX;j++) {
    double now=d.mjdmin+(d.mjdmax-d.mjdmin)*j/NMAX;
    velocity(orb,now,site,&v,&az,&el);
    if(v0<=0 && v>=0 && el>0) {
      double lo=prev,hi=now;
      for(int k=0;k<35;k++) {
        double mid=(lo+hi)/2;velocity(orb,mid,site,&v,&az,&el);
        if(v<0)lo=mid;else hi=mid;
      }
      return (lo+hi)/2;
    }
    prev=now;v0=v;
  }
  return NAN;
}
static int native_plot(const char *device,const char *mode,double rms,double tca) {
  if(cpgopen(device)<=0)return 0;
  cpgscr(0,1,1,1);cpgscr(1,0,0,0);cpgpap(11.,.72);cpgsch(.85);
  char title[160],date[32],epoch[32],xlabel[100];
  mjd2nfd(d.mjdmin,date);mjd2nfd(doy2mjd(orb.ep_year,orb.ep_day),epoch);
  snprintf(xlabel,sizeof(xlabel),"Seconds since %s UTC",date);
  snprintf(title,sizeof(title),"STRF / NORAD %d / %s / RMS %.3f Hz",orb.satno,mode,rms);
  double span=(d.mjdmax-d.mjdmin)*86400,lo=1e30,hi=-1e30,rmax=1;
  for(int i=0;i<d.n;i++) {
    double y=(d.p[i].freq-d.f0),p=(d.p[i].freq0-d.f0);
    lo=fmin(lo,fmin(y,p));hi=fmax(hi,fmax(y,p));rmax=fmax(rmax,fabs(d.p[i].res*1000));
  }
  double pad=fmax((hi-lo)*.08,.01);
  cpgsvp(.10,.95,.43,.86);cpgswin(0,span,lo-pad,hi+pad);
  cpgbox("BCST",0,0,"BCNST",0,0);
  char ylabel[100];snprintf(ylabel,sizeof(ylabel),"Frequency - %.3f MHz (kHz)",d.f0/1000);
  cpglab("",ylabel,"");
  cpgmtxt("T",1.3,0,0,title);
  snprintf(title,sizeof(title),"T_EP %s UTC; RF %.6f MHz; N %d",epoch,d.ffit/1000,d.n);
  cpgmtxt("T",3.1,0,0,title);
  if(isfinite(tca)) {mjd2nfd(tca,date);snprintf(title,sizeof(title),"T_CA %s UTC",date);}
  else snprintf(title,sizeof(title),"T_CA unavailable inside sampled arc");
  cpgmtxt("T",4.9,0,0,title);
  for(int i=0;i<d.n;i++) {
    float x=(d.p[i].mjd-d.mjdmin)*86400;
    cpgsci(4);cpgpt1(x,d.p[i].freq-d.f0,17);
    cpgsci(2);cpgpt1(x,d.p[i].freq0-d.f0,2);
  }
  cpgsci(1);cpgmtxt("B",1,0,0,"Blue: measurements; red crosses: serialized TLE prediction");
  cpgsvp(.10,.95,.12,.32);cpgswin(0,span,-rmax*1.15,rmax*1.15);
  cpgbox("BCNST",0,0,"BCNST",0,0);cpglab(xlabel,"Observed - predicted (Hz)","");
  cpgsls(2);cpgmove(0,0);cpgdraw(span,0);cpgsls(1);
  cpgsci(4);
  for(int i=0;i<d.n;i++)cpgpt1((d.p[i].mjd-d.mjdmin)*86400,d.p[i].res*1000,17);
  cpgsci(1);cpgclos();return 1;
}
int main(int argc,char **argv) {
  if(argc!=7 && argc!=8){fprintf(stderr,"data catalog output residuals mode fit_cut_mjd [plot_device]\n");return 2;}
  char *end;double cut=strtod(argv[6],&end);int ntrain=0;
  if(end==argv[6]||*end||!isfinite(cut)){fprintf(stderr,"Invalid fit cutoff\n");return 2;}
  const char *mode=argv[5];
  int fixed=!strcmp(mode,"seed")||!strcmp(mode,"fixed");
  int ma=!strcmp(mode,"ma")||!strcmp(mode,"mamm")||!strcmp(mode,"maraan")||!strcmp(mode,"mammraan");
  int mm=!strcmp(mode,"mamm")||!strcmp(mode,"mammraan");
  int raan=!strcmp(mode,"raan")||!strcmp(mode,"maraan")||!strcmp(mode,"mammraan");
  if(!fixed&&!ma&&!raan){fprintf(stderr,"Unknown fit mode\n");return 2;}
  d=read_data(argv[1],0,0);d.fitfreq=1;
  for(int i=0;i<d.n;i++){d.p[i].flag=d.p[i].mjd<cut?2:1;if(d.p[i].flag==2)ntrain++;}
  if(ntrain<3){fprintf(stderr,"Insufficient fit data\n");return 2;}
  tle_array_t *catalog=load_tles(argv[2]);
  tle_t *tle=(catalog && catalog->number_of_elements==1)?get_tle_by_index(catalog,0):NULL;
  if(!tle){fprintf(stderr,"Seed catalog must contain exactly one TLE\n");return 2;}
  orb=tle->orbit;d.satname=tle->name;
  double a[7]={orb.eqinc*R2D,orb.ascn*R2D,orb.ecc,orb.argp*R2D,orb.mnan*R2D,orb.rev,orb.bstar};
  double da[7]={0,0,0,0,0,0,0};
  if(ma)da[4]=5.;
  if(mm)da[5]=.1;
  if(raan)da[1]=1.;
  if(!fixed)versafit(ntrain,7,a,da,chisq,0.,1e-5,"n");
  if(!fixed && STRF_FIT_STATUS<0){fprintf(stderr,"Optimizer did not converge\n");return 3;}
  double objective=chisq(a);double fitted_frequency=d.ffit;
  if(!isfinite(objective)||!isfinite(fitted_frequency)||fitted_frequency<=0||orb.rev<=0){fprintf(stderr,"Invalid fit result\n");return 3;}
  char line1[70],line2[70];format_tle(orb,line1,line2);
  FILE *out=fopen(argv[3],"wx");if(!out){perror("Output already exists or cannot be created");return 2;}
  FILE *res=fopen(argv[4],"wx");if(!res){perror("Residual output");fclose(out);remove(argv[3]);return 2;}
  fprintf(out,"NORAD %d\n%s\n%s\n",orb.satno,line1,line2);fclose(out);
  /* Evaluate the serialized TLE, preserving the frequency estimated on selected fit points. */
  tle_array_t *saved=load_tles(argv[3]);tle_t *serialized=saved?get_tle_by_catalog_id(saved,orb.satno):NULL;
  if(!serialized){fclose(res);remove(argv[3]);remove(argv[4]);return 3;}
  orb=serialized->orbit;d.ffit=fitted_frequency;
  if(init_sgdp4(&orb)==SGDP4_ERROR){fclose(res);remove(argv[3]);remove(argv[4]);return 3;}

  fprintf(res,"mjd,site_id,split,observed_hz,predicted_hz,residual_hz,elevation_deg\n");
  double train_ss=0,hold_ss=0,min_alt=90;int nhold=0;
  for(int i=0;i<d.n;i++){
    double v,azi,alt;velocity(orb,d.p[i].mjd,d.p[i].s,&v,&azi,&alt);
    double pred=(1-v/C)*d.ffit,err=(d.p[i].freq-pred)*1000.;
    if(!isfinite(err)||!isfinite(alt)){fclose(res);remove(argv[3]);remove(argv[4]);return 3;}
    d.p[i].freq0=pred;d.p[i].res=err/1000.;
    if(alt<min_alt)min_alt=alt;
    if(d.p[i].flag==2)train_ss+=err*err;else{hold_ss+=err*err;nhold++;}
    fprintf(res,"%.10f,%d,%s,%.6f,%.6f,%.6f,%.6f\n",d.p[i].mjd,d.p[i].site_id,d.p[i].flag==2?"fit":"evaluation",d.p[i].freq*1000,pred*1000,err,alt);
  }
  fclose(res);
  printf("frequency_hz=%.6f rms_hz=%.6f evaluation_rms_hz=%.6f minimum_elevation_deg=%.6f nfit=%d nevaluation=%d\n",d.ffit*1000,sqrt(train_ss/ntrain),nhold?sqrt(hold_ss/nhold):NAN,min_alt,ntrain,nhold);
  char epoch[32],tca_text[32]="unavailable";
  double tca=pass_tca(d.p[0].s);
  mjd2nfd(doy2mjd(orb.ep_year,orb.ep_day),epoch);
  if(isfinite(tca))mjd2nfd(tca,tca_text);
  printf("mode=%s TCA_UTC=%s TEP_UTC=%s tca_mjd=%.10f\n",mode,tca_text,epoch,tca);
  if(argc==8 && !native_plot(argv[7],mode,sqrt((train_ss+hold_ss)/d.n),tca))return 4;
  free(d.p);free_tles(saved);free_tles(catalog);return 0;
}
