#include "header.h"

int main(int argc, char *argv[]){


  bool OVERDENSITY=false, START=false, PEEBLES=false, MASK=false, CALCILMJLM=false, READILMJLM=false, GALAXYCOUNTS=false;
  string inputfile, outputfile, maskfile, IlmJlmfile, mapfile;
  string inputfileroot, outputfileroot;

  int arg=1;
  int lmax=-100;
  int Nmaps=1;
  int nside=512;

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
      case 'P':
	PEEBLES=true;
	arg++;
	break;
      case 'o':
	OVERDENSITY=true;
	arg++;
	break;
      case 'm':
	MASK=true;
	maskfile=string(argv[arg+1]);
	cout<<maskfile<<endl;
	arg+=2;
	break;
      case 'M':
	mapfile=string(argv[arg+1]);
	arg+=2;
	cout<<mapfile<<endl;
	break;
      case 'G':
	GALAXYCOUNTS=true;;
	arg++;
	break;
      case 'R':
	READILMJLM=true;
	IlmJlmfile=string(argv[arg+1]);
	cout<<IlmJlmfile<<endl;
	arg+=2;
	break;
      case 'C':
	CALCILMJLM=true;
	IlmJlmfile=string(argv[arg+1]);
	arg+=2;
	break;
      case 'N':
	nside=atoi(argv[arg+1]);
	arg+=2;
	break;
       case 'L':
	lmax=atoi(argv[arg+1]);
	arg+=2;
	break;
      case 'h':


	cout<<BOLDRED<<"Options:"<<RESET<<endl;
	cout<<endl;
	cout<<YELLOW<<"-I"<<RESET<<"nput: name of the "<<BLUE<<"'inputfile'"<<RESET<<" in ascii"<<endl;
	cout<<YELLOW<<"-O"<<RESET<<"utput: name of the "<<BLUE<<"'outputfile'"<<RESET<<" in ascii format"<<endl;
	cout<<YELLOW<<"-P"<<RESET<<"eebles: Peebles approximation of Cl is calculated; -R or -C and -m as additional option needed!"<<endl;
	cout<<YELLOW<<"-m"<<RESET<<"ask: "<<BLUE<<"'inputfile'"<<RESET<<" of mask in Healpix fits format"<<endl;
	cout<<YELLOW<<"-o"<<RESET<<"verdensity:  overdensity input Alm and Map for Peebles"<<endl;
	cout<<YELLOW<<"-G"<<RESET<<"alaxycounts: if -G the shotnoise of galaxy counts is removed. So input map should be density of counts."<<endl;
	cout<<YELLOW<<"-R"<<RESET<<"eadIlmJlm: "<<BLUE<<"'inputfile'"<<RESET<<" of Ilm and Jlm in ascii format"<<endl;
	cout<<YELLOW<<"-C"<<RESET<<"alcIlmJlm: Ilm and Jlm are calculated and written to "<<BLUE<<"'outputfile'"<<RESET<<";"<<endl
	<<RED<<"   Not recommended! Use CalcIlmJlm before and then use -R option"<<RESET<<endl;
	cout<<YELLOW<<"-M"<<RESET<<"ap: "<<BLUE<<"'mapfile'"<<RESET<<" of map in fits format; only for Peebles approximation needed"<<endl;
	cout<<YELLOW<<"-N"<<RESET<<"side:"<<BLUE<<"'nside'"<<RESET<<" of the map "<<endl;
	cout<<YELLOW<<"-L"<<RESET<<"max:"<<BLUE<<"'lmax'"<<RESET<<" of the output Cl"<<endl;
	cout<<MAGENTA<<"Example:"<<endl<<"./Alm2Cl -I input_alm.fits -O output_cl.dat -m mask.fits -R Ilm_Jlm_file.dat -P -M map.fits -L 100 -N 128"<<endl
	<<" computes the Peebles approximation of the power spectrum of the input alms 'input_alm.fits' with mask 'mask.fits' and writes it to output file 'output_cl.dat';"<<RESET<<endl<<endl;
	cout<<GREEN<<"Queries should be directed to wolz.laura@gmail.com"<<RESET<<endl;
	START=false;
	arg++;
	break;
      default:
	cout<<"bad argument: do Map2Cl -h to get man page. "<<RED<<"Exit now"<<RESET<<endl;
	exit(-1);
	break;
      }
    }else{
      cout<<"bad argument: do Map2Cl -h to get man page. "<<RED<<"Exit now."<<RESET<<endl;
      exit(-1);
    }
  }
  if(!inputfile.empty() && !outputfile.empty()) START=true;

  //constructor for class where analysis is computed
  spherical_shell Shell;
  cout<<"Peebles: "<<PEEBLES<<"   Mask:"<<MASK<<" "<<maskfile<<endl;
  if(START){
   
      //calc/read in of: mask, Ilm, Jlm ans first map
      Shell.set( inputfile, nside, lmax);
      if(MASK){
	cout<<"read mask now!"<<endl;
	Shell.read_mask(maskfile);
      }
      if(READILMJLM){
	cout<<"read Ilm Jlm now!"<<endl;
	Shell.read_I_lm_J_lm(IlmJlmfile);
      }
      if(CALCILMJLM) Shell.calc_I_lm_J_lm(IlmJlmfile);
      
      //read in Alm
      Shell.read_alm(inputfile);     
      if(PEEBLES){
	if(!MASK && (!READILMJLM || !CALCILMJLM)){
	  cout<<"Either mask or IlmJlm are not defined! Cannot calculate Peebles approximation, see  Alm2Cl -h to get man page. "<<RED<<"Exit now."<<RESET<<endl;
	  exit(-1);
	}
	cout<<"Calc PEEbles now!"<<endl;
	Shell.calc_PeeblesCl(OVERDENSITY,  outputfile, mapfile, GALAXYCOUNTS);
      }else{
	Shell.calc_cl(MASK, outputfile, mapfile);
      }
      
    
  }else cout<<"Input or Output file name empty; nothing to do;"<<endl;
    
 return 0;
}
