#include "header.h"

int main(int argc, char *argv[]){


  bool OVERDENSITY=false, START=false, MULTIPLE=false, MASK=false;
  string inputfile, outputfile;
  string inputfileroot, outputfileroot;

  int arg=1;
  int lmax=-100;
  int Nmaps=1;
  double psi=0.0, phi=0.0, theta=0.0;

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
      case 'A':
	psi=atof(argv[arg+1]);
	theta=atof(argv[arg+2]);
	phi=atof(argv[arg+3]);
	arg+=4;
	break;
      case 'h':
	cout<<BOLDRED<<"Options:"<<RESET<<endl;
	cout<<endl;
	cout<<YELLOW<<"-I"<<RESET<<"nput: name of the "<<BLUE<<"'inputfile'"<<RESET<<" in fits"<<endl;
	cout<<YELLOW<<"-O"<<RESET<<"utput: name of the "<<BLUE<<"'outputfile'"<<RESET<<" in fits format"<<endl;
	cout<<YELLOW<<"-A"<<RESET<<"ngles: "<<BLUE<<"'psi' 'theta' 'phi'"<<RESET<<" the Euler angles of the rotation in degrees"<<endl;
	cout<<MAGENTA<<"Example:"<<endl<<"./RotateAlm -I input_alm.fits -O output_cl.dat -A 45.0 45.0 0.0 "<<endl<<" computes the Alm rotated around and writes it to output file 'output_alm.fits';"<<RESET<<endl<<endl;
	cout<<GREEN<<"Queries should be directed to wolz.laura@gmail.com"<<RESET<<endl;
	START=false;
	arg++;
	break;
	
      default:
	cout<<"bad argument: do RotateAlm -h to get man page. "<<RED<<"Exit now"<<RESET<<endl;
	exit(-1);
	break;
      }
    }else{
      cout<<"bad argument: do RotateAlm -h to get man page. "<<RED<<"Exit now."<<RESET<<endl;
      exit(-1);
    }
  }
  if(!inputfile.empty() && !outputfile.empty()) START=true;

  //constructor for class where analysis is computed
  spherical_shell Shell;

  if(START){
      Shell.set(inputfile);
      rotate_alm(Shell.alm, psi*pi/180.0, theta*pi/180.0, phi*pi/180.0);
      Shell.write_alm(outputfile);
     
  }else cout<<"Input or Output file name empty; nothing to do;"<<endl;
    
 return 0;
}
