#include "header.h"

spherical_shell::spherical_shell(){}
spherical_shell::~spherical_shell(){}

const double EPS=1.0e-10; 
const double EPSREL=1.0e-9;
const size_t n=10000;
const double LIGHT=3.0e5;


void spherical_shell::set(int l, string inputfile){

  read_Healpix_map_from_fits(inputfile, this->hp_map, 1, 2);
  this->scheme=hp_map.Scheme();
  this->nside=hp_map.Nside();
  this->order=log2(this->nside);
  this->hp_mask.Set(this->order, this->scheme);
  if(l==-100) this->lmax=2*this->nside+1;
  else this->lmax=l+1;
  this->mmax=this->lmax;
  this->alm.Set(this->lmax, this->mmax);
  this->gal_number=0.0;

  cout<<CYAN<<"ORDER: "<<this->order<<"  NSIDE: "<<this->nside<<"  LMAX: "<<this->lmax<<RESET<<endl; 

}
void spherical_shell::set(int ns, int l){

  this->scheme=RING;
  this->nside=ns;
  this->order=log2(this->nside);
  this->hp_mask.Set(this->order, this->scheme);
  if(l==-100) this->lmax=2*this->nside;
  else this->lmax=l;
  this->mmax=this->lmax;
  this->alm.Set(this->lmax, this->mmax);
  this->gal_number=0.0;
  this->hp_map.Set(this->order, this->scheme);

  cout<<CYAN<<"ORDER: "<<this->order<<"  NSIDE: "<<this->nside<<"  LMAX: "<<this->lmax<<RESET<<endl; 

}

void spherical_shell::set(string inputfile, int ns){

  this->scheme=RING;
  this->nside=ns;
  this->order=log2(this->nside);
  this->hp_map.Set( this->order, this->scheme);
  this->hp_mask.Set(this->order, this->scheme);
  get_almsize(inputfile, this->lmax, this->mmax,2);
  this->alm.Set(this->lmax, this->mmax);
  read_Alm_from_fits(inputfile,this->alm,  this->lmax, this->mmax,2); 
  this->gal_number=0.0;

  cout<<CYAN<<"ORDER: "<<this->order<<"  NSIDE: "<<this->nside<<"  LMAX: "<<this->lmax<<RESET<<endl; 

}

void spherical_shell::set(string inputfile, int ns, int l){

  this->scheme=RING;
  this->nside=ns;
  this->order=log2(this->nside);
  this->hp_map.Set( this->order, this->scheme);
  this->hp_mask.Set(this->order, this->scheme);
  get_almsize(inputfile, this->lmax, this->mmax,2);
  this->alm.Set(this->lmax, this->mmax);
  read_Alm_from_fits(inputfile,this->alm,  this->lmax, this->mmax,2); 
  this->gal_number=0.0;
  this->lmax=l;
  this->mmax=this->lmax;
  cout<<CYAN<<"ORDER: "<<this->order<<"  NSIDE: "<<this->nside<<"  LMAX: "<<this->lmax<<RESET<<endl; 

}





void spherical_shell::set(string inputfile){

  get_almsize(inputfile, this->lmax, this->mmax,2);
  this->alm.Set(this->lmax, this->mmax);
  read_Alm_from_fits(inputfile,this->alm,  this->lmax, this->mmax,2);
  this->mmax=this->lmax;
  this->gal_number=0.0;
  cout<<CYAN<<"  LMAX: "<<this->lmax<<RESET<<endl; 

}

void spherical_shell::read_alm(string inputfile){
  read_Alm_from_fits(inputfile,this->alm,  this->lmax, this->mmax,2);

}
void spherical_shell::read_map(string inputfile){
  read_Healpix_map_from_fits(inputfile, this->hp_map, 1, 2);
}
void spherical_shell::write_map(string outputfile){
  write_Healpix_map_to_fits(outputfile, this->hp_map, PLANCK_FLOAT64);
}

void spherical_shell::write_alm(string outputfile){
  write_Alm_to_fits(outputfile, this->alm, this->lmax,this->mmax, PLANCK_FLOAT64);
}

void spherical_shell::read_mask(string maskfile){

 this->pixelpos=new int[this->hp_map.Npix()];
 this->mask_area=0;
 Healpix_Map<double> buffer;
 read_Healpix_map_from_fits<float64>(maskfile,buffer, 1, 2 );
 this->hp_mask.Import(buffer,false);
 for(int i=0; i<this->hp_mask.Npix(); i++){
    if(this->hp_mask[i]!=0.0){
      this->pixelpos[this->mask_area]=i;
      this->mask_area++;
    }
  }
 cout<<"READ MASK"<<endl;
}




void spherical_shell::calc_alm(bool MASK){

  pointing point_here;
  double ylm;
  if(MASK){
  
      arr<double>weight(2*hp_map.Nside());
    for(int i=0; i<2*hp_map.Nside(); i++) weight[i]=1.0;
    Healpix_Map<double> buffer(this->order, RING);
    buffer.fill(0.0);
    for (int i=0; i<this->hp_map.Npix(); i++) buffer[i]=this->hp_mask[i]*this->hp_map[i];
    map2alm_iter(buffer, this->alm , 10, weight);
   
  }else{
    arr<double>weight(2*hp_map.Nside());
    for(int i=0; i<2*hp_map.Nside(); i++) weight[i]=1.0;
    map2alm_iter(this->hp_map, this->alm ,10, weight);
 
      }
}





//==========================================================================
//==============================Cl==========================================



void spherical_shell::calc_cl(bool MASK, string outputfile, string mapfile){

  arr<double> temp_ps(this->lmax+1);
  double skyfrac;
   read_Healpix_map_from_fits(mapfile, this->hp_map, 1,2);
  PowSpec powspec(1, this->lmax);
  if(MASK){
    skyfrac=double(this->mask_area)/double(this->hp_map.Npix());
    
  }else{
    skyfrac=1.0;
  }
  for(int i=0; i<this->hp_map.Npix();i++) this->gal_number+=this->hp_map[i];

  double shotnoise=1.0/(this->gal_number/(4.0*pi*skyfrac));
  cout<<"Skyfraction: "<<skyfrac<<endl;
  /* arr<double>winpix(4*this->nside);
  string winpix_directory="/home/lwolz/Libraries/Healpix_donatello/Healpix_2.20a/data/";
  read_pixwin(winpix_directory,this->nside, winpix);
  */

  cout<<"Shotnoise: "<<shotnoise<<endl;

  ofstream file(&outputfile[0]);
  if(file.is_open()){
    file<<"RMS:    Average: "<<endl;
    file<<hp_map.rms()<<"   "<<hp_map.average()<<endl;
    for(int l=0; l<this->lmax; l++){
      temp_ps[l]=0.0;
      for(int m=0; m<=l; m++){
	temp_ps[l]+=( pow(this->alm(l,m).re, 2.0) + pow(this->alm(l,m).im, 2.0) ); 
	if(m==0) temp_ps[l]/=2.0;   
      }
       
      temp_ps[l]/=((double(l)+0.5)*skyfrac);
      // temp_ps[l]/=(winpix[l]*winpix[l]);
      //temp_ps[l]/=pow(this->hp_map.average(),2.0) ;
      //temp_ps[l]-=shotnoise;
      file<<l<<setw(15)<<temp_ps[l]<<endl;
    }
  }else{
    cout<<"Cannot open outputfile! "<<RED<<" Exit now."<<RESET<<endl;
    exit(-1);
  }
  this->powspec.Set(temp_ps);
}





void spherical_shell::calc_PeeblesCl(bool OVERDENSITY, string outputfile, string mapfile, bool GALAXYCOUNTS){

  arr<double> Cltemp(this->lmax+1);

  double deltaomega=0.0;
  double complexabs=0.0;

  read_Healpix_map_from_fits(mapfile, this->hp_map, 1,2);
//  arr<double>winpix(4*this->nside);
//  string winpix_directory="/home/lwolz/Library/Healpix_2.20a/data/";
//  read_pixwin(winpix_directory,this->nside, winpix);

  deltaomega=double(this->mask_area)*4.0*PI/double(this->hp_map.Npix());
  this->gal_number=0.0;
  
  for(int i=0; i<this->hp_map.Npix(); i++){
    if(this->hp_mask[i]==1) this->gal_number+=this->hp_map[i];
  }
  double  skyfrac=double(this->mask_area)/double(this->hp_map.Npix());
  double shotnoise=1.0/(this->gal_number/(4.0*pi*skyfrac));
  double omega_pix=4.0*pi/this->hp_map.Npix();
  double factor=this->gal_number/deltaomega*omega_pix;


  cout<<"galaxy N: "<<this->gal_number<<endl<<"mask area: "<<this->mask_area<<endl;
  if(GALAXYCOUNTS) cout<<"Shotnoise: "<<shotnoise<<endl;

  double testsum=0.0;
  for(int i=0; i<this->lmax; i++) Cltemp[i]=0.0;
  double average=this->gal_number/this->mask_area;
  ofstream file(&outputfile[0]);
  if(file.is_open()){
    file<<"RMS:    Average: "<<endl;
    file<<hp_map.rms()<<"   "<<average<<endl;

    if(OVERDENSITY){
      for(int l=0; l<this->lmax; l++){
	for(int m=0; m<=l; m++){
	  Cltemp[l]+=( (pow(this->alm(l,m).re,2.0)+pow(this->alm(l,m).im, 2.0))/ this->J[l][m]);
	  if(m==0) Cltemp[l]/=2.0;
	}
	Cltemp[l]/=(double(l)+0.5) ;
	file<<l<<setw(15)<<Cltemp[l]<<endl;
      }
    }else{
      for(int l=0; l<this->lmax; l++){
	for(int m=0; m<=l; m++){
	  complexabs=pow(this->alm(l,m).re- factor*this->I[l][m].re, 2.0)+pow(this->alm(l,m).im- factor*this->I[l][m].im, 2.0);
	  Cltemp[l]+= (complexabs/ this->J[l][m]);
      	  if(m==0) Cltemp[l]/=2.0;

	}
	Cltemp[l]/=(double(l)+0.5);
	//	Cltemp[l]/=pow(average,2.0);
	if(GALAXYCOUNTS) Cltemp[l]-=shotnoise;
	//Cltemp[l]/=(winpix[l]*winpix[l]);

	file<<l<<setw(15)<<Cltemp[l]<<endl;
      }
    }
  }else{
    cout<<"Cannot open outputfile! "<<RED<<" Exit now."<<RESET<<endl;
    exit(-1);
  }
  this->Peebles.Set(Cltemp);

}




//============================================================================
//================   Ilm Jlm ===============================================

void spherical_shell::calc_I_lm_J_lm(string IlmJlmfile){

  this->I=new xcomplex<double>*[this->lmax];
  this->J=new double*[this->lmax];
  for(int i=0; i<this->lmax; i++){
    this->I[i]=new xcomplex<double>[this->lmax];
    this->J[i]=new double[this->lmax];
  }
  pointing pointer;
  double ylm=0.0;
  xcomplex<double> Ylm;
  double pixarea=4.0*PI/double(this->hp_map.Npix());
  double deltaomega=double(this->mask_area)*4.0*PI/double(this->hp_map.Npix());
  for(int l=0; l<this->lmax; l++){
      for(int m=0; m<=l; m++){
	this->I[l][m].re=0.0;
	this->I[l][m].im=0.0;
	this->J[l][m]=0.0;

      }
  }

  cout<<IlmJlmfile<<endl;
  ofstream file(&IlmJlmfile[0]);
  if(file.is_open()){
for(int l=0; l<this->lmax; l++){
  cout<<l<<"   ";
      for(int m=0; m<=l; m++){
	for(int i=0; i<this->mask_area; i++){
	  pointer=this->hp_map.pix2ang(double(this->pixelpos[i]));
	  ylm=gsl_sf_legendre_sphPlm(double(l), abs(double(m)) ,cos(pointer.theta));
	  Ylm.re=ylm*cos(double(m)*pointer.phi);
	  Ylm.im=ylm*sin(double(m)*pointer.phi);     
	  //if(m<0) Ylm=pow(-1.0, double(m))*Ylm.conj();
	  this->I[l][m].re+=(Ylm.re*pixarea);
	  this->I[l][m].im+=(-Ylm.im*pixarea);	
	  this->J[l][m]+=(pow(abs(Ylm), 2.0)*pixarea);  
	}

      }
 }
 
    for(int l=0; l<this->lmax; l++){
      for(int m=0; m<=l; m++){
	file<<l<<"  "<<m<<"  "<<I[l][m].re<<"   "<<I[l][m].im<<"   "<<J[l][m]<<endl;
      }
    }
    file.close();
  }else{
    cout<<"Cannot write to file in Ilm and Jlm! "<<RED<<" Exit now."<<RESET<<endl;
    exit(-1);
  }
}





void spherical_shell::calc_I_lm_J_lm(string IlmJlmfile, int startl, int stopl){

  this->I=new xcomplex<double>*[this->lmax];
  this->J=new double*[this->lmax];
  for(int i=0; i<this->lmax; i++){
    this->I[i]=new xcomplex<double>[this->lmax];
    this->J[i]=new double[this->lmax];
  }
  pointing pointer;
  double ylm=0.0;
  xcomplex<double> Ylm;
  double pixarea=4.0*PI/double(this->hp_map.Npix());
  double deltaomega=double(this->mask_area)*4.0*PI/double(this->hp_map.Npix());
  for(int l=0; l<this->lmax; l++){
      for(int m=0; m<=l; m++){
	this->I[l][m].re=0.0;
	this->I[l][m].im=0.0;
	this->J[l][m]=0.0;

      }
  }

  cout<<IlmJlmfile<<endl;
  ofstream file(&IlmJlmfile[0]);
  if(file.is_open()){
for(int l=startl; l< stopl; l++){
  cout<<l<<"   ";
      for(int m=0; m<=l; m++){
	for(int i=0; i<this->mask_area; i++){
	  pointer=this->hp_map.pix2ang(double(this->pixelpos[i]));
	  ylm=gsl_sf_legendre_sphPlm(double(l), abs(double(m)) ,cos(pointer.theta));
	  Ylm.re=ylm*cos(double(m)*pointer.phi);
	  Ylm.im=ylm*sin(double(m)*pointer.phi);     
	  //if(m<0) Ylm=pow(-1.0, double(m))*Ylm.conj();
	  this->I[l][m].re+=(Ylm.re*pixarea);
	  this->I[l][m].im+=(-Ylm.im*pixarea);	
	  this->J[l][m]+=(pow(abs(Ylm), 2.0)*pixarea);  
	}

      }
 }
 
    for(int l=startl; l<stopl ; l++){
      for(int m=0; m<=l; m++){
	file<<l<<"  "<<m<<"  "<<this->I[l][m].re<<"   "<<this->I[l][m].im<<"   "<<this->J[l][m]<<endl;
      }
    }
    file.close();
  }else{
    cout<<"Cannot write to file in Ilm and Jlm! "<<RED<<" Exit now."<<RESET<<endl;
    exit(-1);
  }
}




void spherical_shell::read_I_lm_J_lm(string IlmJlmfile){
  
  ifstream file(&IlmJlmfile[0]);
  int templ, tempm;
  this->I=new xcomplex<double>*[this->lmax];
  this->J=new double*[this->lmax];
  for(int i=0; i<this->lmax; i++){
    this->I[i]=new xcomplex<double>[this->lmax];
    this->J[i]=new double[this->lmax];
  }
  if(file.is_open()){
    for(int l=0; l<this->lmax; l++){
      for(int m=0; m<=l; m++){
	file>>templ>>tempm>>I[templ][tempm].re>>I[templ][tempm].im>>J[templ][tempm];
      }
    }
  }else{
    cout<<"Cannot read Ilm and Jlm! "<<RED<<" Exit now."<<RESET<<endl;
    exit(-1);
  }
}


//-----------------------------------------------------------------------------------------------------------------------------------------------------
//CREATE DIFFERENT KIND OF MASKS
void spherical_shell::calc_mask_lognorm(string outputfile){
 
  pointing pointer;
  this->mask_area=0;
  for(int i=0; i<this->hp_mask.Npix(); i++){
    pointer=this->hp_mask.pix2ang(i);
    if(( pointer.phi<=(270.0*PI/180.0) &&  (pointer.phi)>=(110.0*PI/180.0) && pointer.theta>=((20.0)*PI/180.0) &&  pointer.theta<=((95.0)*PI/180.0))){
      this->hp_mask[i]=1.0;
      this->mask_area++;
    }else{
      this->hp_mask[i]=0.0;      
    }
  }
write_Healpix_map_to_fits(outputfile, this->hp_mask, PLANCK_FLOAT64);
}


void spherical_shell::calc_mask_halfsky(string outputfile){
 
  pointing pointer;
  this->mask_area=0;
  for(int i=0; i<this->hp_mask.Npix(); i++){
    pointer=this->hp_mask.pix2ang(i);
    if(( pointer.theta<=pi/2)  ){//|| (2.0*PI-pointer.phi)<=this->phi_lim )){
      this->hp_mask[i]=1.0;
      this->mask_area++;
    }else{
      this->hp_mask[i]=0.0;      
    }
  }
write_Healpix_map_to_fits(outputfile, this->hp_mask, PLANCK_FLOAT64);
}


void spherical_shell::calc_mask_fg(string outputfile, string fgfile, double tlim){
  Healpix_Map <double> fg;
  read_Healpix_map_from_fits(fgfile, fg, 1,2);
  this->mask_area=0;
  this->hp_mask.fill(0.0);
  for(int i=0; i<fg.Npix(); i++){
    if(fg[i]<=tlim) this->hp_mask[i]=1.0;
  }
write_Healpix_map_to_fits(outputfile, this->hp_mask, PLANCK_FLOAT64);
}




void spherical_shell::calc_mask(string outputfile, double *theta, double *phi){
 
  pointing pointer;
  this->mask_area=0;
  for(int i=0; i<this->hp_mask.Npix(); i++){
    pointer=this->hp_mask.pix2ang(i);
    if(( pointer.theta>=theta[0] && pointer.theta<=theta[1]) && (pointer.phi>=phi[0] && pointer.phi<=phi[1])  ) {
      this->hp_mask[i]=1.0;
      this->mask_area++;
    }else{
      this->hp_mask[i]=0.0;      
    }
  }
  write_Healpix_map_to_fits<float64>(outputfile, this->hp_mask, PLANCK_FLOAT64);

}







void spherical_shell::winpix( int whichCl){
 
  arr<double>winpix(4*this->nside);
  string winpix_directory="/home/lwolz/Libraries/Healpix_leo/Healpix_2.20a/data/";
  read_pixwin(winpix_directory,this->nside, winpix);
  //for(int i=0; i<4*this->nside; i++) cout<<winpix[i]<<"  "<<pow(winpix[i],2.0)<<endl;
  for(int l=0; l<this->lmax; l++){
    if(whichCl==0) this->powspec.tt(l)/=(winpix[l]*winpix[l]);
    if(whichCl==1) this->Peebles.tt(l)/=(winpix[l]*winpix[l]);
    if(whichCl==2) this->mix_Cl.tt(l)/=(winpix[l]*winpix[l]);
  } 

}
