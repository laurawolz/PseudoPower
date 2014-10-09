#include "header.h"

int main(int argc, char *argv[]){

bool OVERDENSITY=false, START=false,  MASK=false;
  string  outputfile, maskfile, IlmJlmfile;

  int arg=1;

  int lmin=0;
  int lmax=-100;
  int nside=512;
 while (arg < argc) {
    if (argv[arg][0] == '-') {
      switch (argv[arg][1]) {
      case 'l':
	lmin=atoi(argv[arg+1]);
	arg+=2;
	break;
      case 'L':
	lmax=atoi(argv[arg+1]);
	arg+=2;
	break;
      case 'O':
 	IlmJlmfile=string(argv[arg+1]);
	arg+=2;
	break;
      case 'm':
	MASK=true;
	maskfile=string(argv[arg+1]);
	cout<<maskfile<<endl;
	arg+=2;
	break;
      case 'N':
	nside=atoi(argv[arg+1]);
	arg+=2;
	break;
      case 'h':
	cout<<RED<<"This function can have long running times for high l and high Nside. It is advised to run it for small l ranges in parallel and copy all output files together using the 'cat' command."<<endl;
	cout<<BOLDRED<<"Options:"<<RESET<<endl;
	cout<<endl;

	cout<<YELLOW<<"-O"<<RESET<<"utput: name of the "<<BLUE<<"'outputfile'"<<RESET<<" in ascii format"<<endl;

	cout<<YELLOW<<"-m"<<RESET<<"ask: "<<BLUE<<"'inputfile'"<<RESET<<" of mask in fits format"<<endl;

	cout<<YELLOW<<"-N"<<RESET<<"side:"<<BLUE<<"'nside'"<<RESET<<" of the mask and map"<<endl;
	cout<<YELLOW<<"-l"<<RESET<<"min: "<<BLUE<<"'lmin'; "<<RESET<<"start with 0 for complete runs;"<<endl;
	cout<<YELLOW<<"-L"<<RESET<<"max: "<<BLUE<<"'lmax'"<<RESET<<endl;
	cout<<MAGENTA<<"Example:"<<endl<<".IlmJlm  -O ilmjlm_mask1.dat -m mask1.fits -N 128 -l 0 -L 100"<<endl<<"; Calculates the Ilm and Jlm for mask1 with resolution Nside=128 for the multipoles 0 to 100"<<RESET<<endl<<endl;

	START=false;
	arg++;
	break;
      default:
	cout<<"bad argument: do IlmJlm -h to get man page. "<<RED<<"Exit now"<<RESET<<endl;
	exit(-1);
	break;
      }
    }else{
      cout<<"bad argument: do IlmJlm -h to get man page. "<<RED<<"Exit now."<<RESET<<endl;
      exit(-1);
    }
  }
  if( !IlmJlmfile.empty()) START=true;

  //constructor for class where analysis is computed
  spherical_shell Shell;
  if(START){
        
    Shell.set(nside,lmax);   
    Shell.read_mask(maskfile);
    cout<<"here"<<endl;
    Shell.calc_I_lm_J_lm(IlmJlmfile,lmin,lmax);
 
 
  }else cout<<" Output file name empty; nothing to do;"<<endl;
    
 return 0;
}
