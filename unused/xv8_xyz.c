#include <stdio.h>
#include <stdlib.h>

void write_xv8_(char*,int*,float*,float*,float*,int);

int main(int argc, char *argv[]) {
  int j, n, PART, dumpi, k;
  double time8, CHARGE, MASS;
  char filename[100], line[300];
  FILE *chargefile, *out;

  if(argc != 4 && argc != 6) {
    fprintf(stderr,"Syntax Error: %s <configuration> <output file> <ion.dat file or # of particles> <charge> <mass>\n",argv[0]);
    return 0;
  }

  if(argc == 6) {
    PART=atoi(argv[3]);
    sscanf(argv[4],"%lf",&CHARGE);
    sscanf(argv[5],"%lf",&MASS);
  }
  else {
    chargefile=fopen(argv[3],"r");
    fgets(line,300,chargefile);
    sscanf(line,"%d",&PART);
  }

  double pos[PART][3], vel[PART][3], charge[PART], mass[PART];

  if(argc == 4) {
    for(j=0;j<PART;j++) {
      fgets(line,300,chargefile);
      sscanf(line,"%d %lf %lf",&dumpi,&charge[j],&mass[j]);
    }
    fclose(chargefile);
  }
  else
    for(j=0;j<PART;j++) {
      charge[j]=CHARGE;
      mass[j]=MASS;
    }

  
  k=strlen(argv[1]);
  read_xv8_(argv[1],&PART,&time8,&pos[0],&vel[0],k);

  out=fopen(argv[2],"w");
  fprintf(out,"%d\n[name]\n",PART);  

  for(j=0;j<PART;j++) {
    if (charge[j]<1.5)
      fprintf(out,"H %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<2.5)
      fprintf(out,"He %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<3.5)
      fprintf(out,"Li %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<4.5)
      fprintf(out,"Be %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<5.5)
      fprintf(out,"B %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<6.5)
      fprintf(out,"C %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<7.5)
      fprintf(out,"N %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<8.5)
      fprintf(out,"O %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<9.5)
      fprintf(out,"F %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<10.5)
      fprintf(out,"Ne %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<11.5)
      fprintf(out,"Na %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<12.5)
      fprintf(out,"Mg %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<13.5)
      fprintf(out,"Al %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<14.5)
      fprintf(out,"Si %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<15.5)
      fprintf(out,"P %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<16.5)
      fprintf(out,"S %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<17.5)
      fprintf(out,"Cl %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<18.5)
      fprintf(out,"Ar %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<19.5)
      fprintf(out,"K %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<20.5)
      fprintf(out,"Ca %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<21.5)
      fprintf(out,"Sc %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<22.5)
      fprintf(out,"Ti %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<23.5)
      fprintf(out,"V %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<24.5)
      fprintf(out,"Cr %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<25.5)
      fprintf(out,"Mn %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<26.5)
      fprintf(out,"Fe %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<27.5)
      fprintf(out,"Co %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<28.5)
      fprintf(out,"Ni %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<29.5)
      fprintf(out,"Cu %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<30.5)
      fprintf(out,"Zn %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<31.5)
      fprintf(out,"Ga %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<32.5)
      fprintf(out,"Ge %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<33.5)
      fprintf(out,"As %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<34.5)
      fprintf(out,"Se %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<35.5)
      fprintf(out,"Br %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<36.5)
      fprintf(out,"Kr %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<37.5)
      fprintf(out,"Rb %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<38.5)
      fprintf(out,"Sr %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<39.5)
      fprintf(out,"Y %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<40.5)
      fprintf(out,"Zr %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<41.5)
      fprintf(out,"Nb %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<42.5)
      fprintf(out,"Mo %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<43.5)
      fprintf(out,"Tc %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<44.5)
      fprintf(out,"Ru %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<45.5)
      fprintf(out,"Rh %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<46.5)
      fprintf(out,"Pd %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<47.5)
      fprintf(out,"Ag %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<48.5)
      fprintf(out,"Cd %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else if(charge[j]<49.5)
      fprintf(out,"In %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
    else
      fprintf(out,"Sn %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
  }

  fclose(out);

}
