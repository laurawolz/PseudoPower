#include "header.h"

int main(int argc, char *argv[]){


  bool OVERDENSITY=false, START=false, MASK=false;
  string inputfile, outputfile, maskfile;
  string inputfileroot, outputfileroot;

  int arg=1;
  int lmax=-100;
  int Nmaps=1;

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
      case 'L':
 	lmax=atoi(argv[arg+1]);
	arg+=2;
	break;
      case 'm':
	MASK=true;
	maskfile=string(argv[arg+1]);
	arg+=2;
	break;
      case 'h':
	cout<<BOLDRED<<"Options:"<<RESET<<endl;
	cout<<endl;
	cout<<YELLOW<<"-I"<<RESET<<"nput: name of the "<<BLUE<<"'inputfile'"<<RESET<<" in fits"<<endl;
	cout<<YELLOW<<"-O"<<RESET<<"utput: name of the "<<BLUE<<"'outputfile'"<<RESET<<" in Healpix fits format"<<endl;
	cout<<YELLOW<<"-L"<<RESET<<"max: "<<BLUE<<"'lmax'"<<RESET<<" maximum multipole of Cl; default lmax=2*Nside of input map"<<endl;
	cout<<YELLOW<<"-m"<<RESET<<"ask: "<<BLUE<<"'inputfile'"<<RESET<<" of mask in fits format; cuts and weights survey window according to mask; "<<endl;
	cout<<MAGENTA<<"Example:"<<endl<<"./Map2Alm -I input_map.fits -O output_alm.fits -L 100 "<<endl<<" computes the Alms of the input map 'input_map.fits'to lmax=100 and writes it to output file 'output_alm.fits';"<<RESET<<endl<<endl;
	cout<<GREEN<<"Queries should be directed to wolz.laura@gmail.com"<<RESET<<endl;

	START=false;
	arg++;
	break;
      default:
	cout<<"bad argument: do Map2Alm -h to get man page. "<<RED<<"Exit now"<<RESET<<endl;
	exit(-1);
	break;
      }
    }else{
      cout<<"bad argument: do Map2Alm -h to get man page. "<<RED<<"Exit now."<<RESET<<endl;
      exit(-1);
    }
  }
  if(!inputfile.empty() && !outputfile.empty()) START=true;

  //constructor for class where analysis is computed
  spherical_shell Shell;

  if(START){
 
      //for first file set parameter and set first map
	Shell.set(lmax, inputfile);
	if(MASK) Shell.read_mask(maskfile);
      
      //calculation of alm
      Shell.calc_alm(MASK);
      Shell.write_alm(outputfile);
    
  }else cout<<"Input or Output file name empty; nothing to do;"<<endl;
    
 return 0;
}
