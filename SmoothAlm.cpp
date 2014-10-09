#include "header.h"

int main(int argc, char *argv[]){


  bool OVERDENSITY=false, START=false,  BANDLIMIT=false;
  string inputfile, outputfile;
  string inputfileroot, outputfileroot;

  int arg=1;
  int lmax=-100;
  double fwhm=0.0;
  int llim=-100;

 while (arg < argc) {
    if (argv[arg][0] == '-') {
      switch (argv[arg][1]) {
      case 'I':
	inputfile=string(argv[arg+1]);
	arg+=2;
	break;
      case 'O':
 	outputfile=string(argv[arg+1]);
	arg+=2;
	break;
      case 'F':
	fwhm=atof(argv[arg+1]);
	arg+=2;
	break;
      case 'B':
	BANDLIMIT=true;
	llim=atof(argv[arg+1]);
	arg+=2;
	break;
      case 'h':
	cout<<BOLDRED<<"Options:"<<RESET<<endl;
	cout<<endl;
	cout<<YELLOW<<"-I"<<RESET<<"nput: name of the "<<BLUE<<"'inputfile'"<<RESET<<" in ascii"<<endl;
	cout<<YELLOW<<"-O"<<RESET<<"utput: name of the "<<BLUE<<"'outputfile'"<<RESET<<" in ascii format"<<endl;
	cout<<YELLOW<<"-F"<<RESET<<"WHM: "<<BLUE<<"'fwhm'"<<RESET<<" of the Gaussian beam to smooth with in degree"<<endl;
	cout<<YELLOW<<"-B"<<RESET<<"andlimit: "<<BLUE<<"'llim'"<<RESET<<", the l after which all Alm are set to zero "<<endl;
	cout<<MAGENTA<<"Example:"<<endl<<"./SmoothAlm -I input_alm.fits -O output_cl.dat -F 1.0 "<<endl<<" computes the Alm smoothed with a Gaussian beam with full width half maximum 1.0 degree and writes it to output file 'output_alm.dat';"<<RESET<<endl<<endl;
	cout<<GREEN<<"Queries should be directed to wolz.laura@gmail.com"<<RESET<<endl;
	START=false;
	arg++;
	break;
      default:
	cout<<"bad argument: do Smoothalm -h to get man page. "<<RED<<"Exit now"<<RESET<<endl;
	exit(-1);
	break;
      }
    }else{
      cout<<"bad argument: do SmoothAlm -h to get man page. "<<RED<<"Exit now."<<RESET<<endl;
      exit(-1);
    }
  }
  if(!inputfile.empty() && !outputfile.empty()) START=true;
 
  //constructor for class where analysis is computed
  spherical_shell Shell;
  cout<<START<<endl;
  if(START){
   Shell.set(inputfile);
      smoothWithGauss(Shell.alm, fwhm*pi/180.0);
      if(BANDLIMIT){
	cout<<" Setting Alm above l="<<llim<<" to zero;"<<endl;
	for (int l=llim;l<Shell.alm.Lmax();l++){
	  for (int m=0;m<=l;m++){
	    Shell.alm(l,m).re=0.0;
	    Shell.alm(l,m).im=0.0;
	  }
	}
      }
      Shell.write_alm(outputfile);
       
  }else cout<<"Input or Output file name empty; nothing to do;"<<endl;
    
 return 0;
}
