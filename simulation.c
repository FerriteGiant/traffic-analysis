//
//
//gcc simulation.c -lfftw3 -lm -o sim.exe
//
//
//#include <complex.h>
#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define PI 3.1415926535897932384626433832795

void distProbDistributionHist(int *posArray, int *velArray, int *distArray);
void asciiDisplay(int *posArray, int *velArray);
void trackDiffusion(int *posArray, int *velArray);
void pairCorrelation(int *posArray, int *velArray, int *distArray);
//void fftSetup(double **track2D,fftw_complex **outHalf,fftw_complex **outFull, fftw_plan plan_fwd, int nx, int ny);
void trackHistory(int *posArray,double *track2D,int counter);
void fft(double *track2D,fftw_complex *outHalf,double *outFull,fftw_plan plan_fwd, int nx, int ny);

void runSimulation(int *posArray, int *velArray);
void updateCar(int *pos1, int *vel1, int *pos2);
float welch(int timeStep);
float Hamm(int timeStep);
void delay(int milliseconds);

//DEFINE PARAMETERS
int numCars;
float density;
int maxVel;
int trackLength = 8192;
int sampleRate = 1;
int dataStartStep = 100000;
int numSteps,numSamples;
int fftSamples = 2048; //Time steps per FFT
int fftRuns; //Number of matrices to do ffts of
char *fileName;
FILE *outputFile;


int main(int argc, char *argv[])
{


if (argc != 5)
{
  printf("Usage: %s <fftRuns> <density> <maxVel> <save directory>\n",argv[0]);
  return(-1);
}
else
{ //First argument
//printf("test");
  fftRuns = atoi(argv[1]);
  numSamples = fftRuns*fftSamples;
  numSteps = (numSamples*sampleRate)+dataStartStep;
  
  //Second argumnet
  density = atof(argv[2]);
  numCars = trackLength*density;

  //Third argument
  maxVel = atoi(argv[3]);

  //Create filename to save to
  //float dataPoints = numSamples*numCars;
  //printf("dataPoints: %f",dataPoints);
  asprintf(&fileName, "%svel%d_density%.2f_fft%d_track%d_runs%.0e"
                          "_numCars%d_halfdata.csv",\
                          argv[4],maxVel,density,fftSamples,trackLength,\
                          (float)fftRuns,numCars);
  
  outputFile = fopen(fileName,"w");

  if (outputFile == 0)
  {
    printf("Could not open file %s.\n",argv[2]);
    return(-1);
  }

  if (!(numSteps > sampleRate*2))
  {
    printf("Number of steps not a valid entry," 
    "or not more than twice the sample rate of %d.\n",sampleRate);
    return(-1);
  }
}

//Check that the number of samples for each fft is equal to or less than
//the total number of samples to be taken. 
if(fftSamples > numSamples)
{printf("fftSamples must be less than or equal to numSamples\n");return(-1);}
//Also check that numSamples is an integer multiple of fftSamples
if(!(numSamples%fftSamples==0))
{printf("numSamples must be an integer multiple of fftSamples\n");return(-1);}
//And check that fftSamples is an integer multiple of sampleRate
if(!(fftSamples%sampleRate==0))
{printf("fftSamples must be an integer multiple of sampleRate\n");return(-1);}


//ALLOCATE DATA ARRAYS
int i,j,k,fftRun,step;
int *posArray, *velArray;
posArray = malloc((numCars+1)*sizeof(int));
velArray = malloc((numCars+1)*sizeof(int));

srand (time(NULL));

//Initialize starting positions and speeds 
for (i=0;i<numCars;i++){
  posArray[i] = i;
  velArray[i] = 0;
  }

// setup stuff
int counter = 0;
int *distNextArray,*distAnyArray;
distNextArray = malloc((trackLength)*sizeof(int));
distAnyArray = malloc((trackLength)*sizeof(int));
for(i=0;i<trackLength;i++){
  distNextArray[i]=0;
  distAnyArray[i]=0;
}

//Setup DFT stuff
double *track2D;
double *outFull;
fftw_complex *outHalf;
fftw_plan plan_fwd;//,plan_bwd; //Note that this creates 'opaque pointers'
int ny = trackLength;
int nx = fftSamples*2;//pad with zeros for reasons
int nyh = ny/2+1;
track2D = fftw_malloc(sizeof(double)*nx*ny);
outHalf = fftw_malloc(sizeof(fftw_complex)*nx*nyh);
outFull = fftw_malloc(sizeof(double)*nx*ny);

for(i=0;i<nx*ny;i++){
  outFull[i]=0;
}
  plan_fwd = fftw_plan_dft_r2c_2d(nx,ny,track2D,outHalf,FFTW_ESTIMATE);
//fftSetup(&track2D,&outHalf,&outFull,plan_fwd,nx,ny);

//  int *pairCorrArray;
//  pairCorrArray = malloc((trackLength)*sizeof(double));
//  track = malloc((trackLength)*sizeof(int));
//  for(i=0;i<trackLength;i++)
//    pairCorrArray[i]=0;

//Initialization runs
for (step=0;step<dataStartStep;step+=sampleRate){
  runSimulation(posArray, velArray);//Will run sampleRate times
}


//Data taking runs
//FFT will happen each time this runs and the results will be averaged
for (fftRun=0;fftRun<fftRuns;fftRun++){
  counter=0;

  //zero out array for holding track history
  for(j=0;j<nx*ny;j++)
    track2D[j]=0;

  for (i=0;i<fftSamples;i++){
    runSimulation(posArray, velArray);//Will run sampleRate times
  
  
    //// CALLS TO VARIOUS ANALYSIS FUNCTIONS TO RUN AT EACH SAMPLE STEP ////
    //distProbDistributionHist(posArray, velArray, distNextArray);
    //asciiDisplay(posArray, velArray);
    //trackDiffusion(posArray, velArray);
    //pairCorrelation(posArray, velArray,distAnyArray);
    trackHistory(posArray,track2D,counter); //Store output of a given step
  
    counter++;
  }
  //Check that the counter has counted up to fftSamples
  if(counter != fftSamples){
    printf("Counter(%d) & fftSamples(%d) mismatch\n",counter,fftSamples);return(-1);
    }

  
  fft(track2D,outHalf,outFull,plan_fwd,nx,ny);

  if((fftRun+1)%100 == 0)
        printf("Completed %d of %d ffts\n",fftRun+1,fftRuns);
}
int index;
double avg;
  //Print result of DFT
  fprintf(outputFile,"numCars,%d,maxVel,%d,trackLength,%d,sampleRate,%d,"
                    "dataStartStep,%d,fftSamples,%d,fftRuns,%d\n",\
                    numCars,maxVel,trackLength,sampleRate,\
                    dataStartStep,fftSamples,fftRuns);
  int maxny = ny/2+1; //Only store the nonredundant part of the data
  for (i=0;i<nx;i++){
    for (j=0;j<maxny;j++){
      index = i*ny+j;
      avg = outFull[index]/((double)fftRuns*(double)fftSamples*(double)trackLength);
      if(j<maxny-1)
        fprintf(outputFile,"%.6f,",avg);
      else
        fprintf(outputFile,"%.6f",avg);
    }
    fprintf(outputFile,"\n");
  }

// PRINT DISTANCE HISTOGRAM
//fprintf(outputFile,"G(r),G_binCount,P(r),P_binCount\n");
//double pofr,gofr; //P(r), G(r)
//for(i=0;i<120;i++){
//  pofr = (float)distNextArray[i]/(counter*numCars);
//  gofr = (float)distAnyArray[i]/(counter*numCars);
//  fprintf(outputFile,"%.8f,%d,%.8f,%d\n",\
//    gofr,distAnyArray[i],pofr,distNextArray[i]);
//}

//fclose(outputFile);
//fftw_destroy_plan(plan_bwd);
fftw_destroy_plan(plan_fwd);
fftw_free(outHalf);
free(outFull);
free(distNextArray);
free(distAnyArray);
fftw_free(track2D);
free(posArray);
free(velArray);
free(fileName);
}/////// END MAIN


/////////////////////////////////////////////////////////////
// FFT
/////////////////////////////////////////////////////////////

//void fftSetup(double **track2D,fftw_complex **outHalf,fftw_complex **outFull, fftw_plan plan_fwd, int nx, int ny){
//  int i;
//  int nyh = ny/2+1;
//  *track2D = fftw_malloc(sizeof(double)*nx*ny);
//  *outHalf = fftw_malloc(sizeof(fftw_complex)*nx*nyh);
//  *outFull = fftw_malloc(sizeof(fftw_complex)*nx*ny);
//  for(i=0;i<nx*ny;i++){
//    **outFull[i][0]=0;
//    **outFull[i][1]=0;
//  }
//
//  plan_fwd = fftw_plan_dft_r2c_2d(nx,ny,*track2D,*outHalf,FFTW_ESTIMATE);
//  //plan_bwd = fftw_plan_dft_c2r_2d(nx,ny,out,in,FFTW_ESTIMATE);
//}

//Fill in array which stores the track state at each sample step
void trackHistory(int *posArray,double *track2D, int counter)
{
  int i,t;
//  double sVal,wavelength,period;
//  wavelength = (double)trackLength/1.0;
//  period = (double)fftSamples/1.0;
  for(i=0;i<numCars;i++){
    //track2D[posArray[i]+counter*trackLength]=welch(counter);
    track2D[posArray[i]+counter*trackLength]=Hamm(counter);
    //printf("%.3f ",welch(counter));
    //printf("%.3f ",Hamm(counter));
  //double sVal,wavelength,period;
  //wavelength = (double)trackLength/1.0;
  //period = (double)fftSamples/1.0;
  }
  //printf("\n");
//    for(i=0;i<trackLength;i++){
//      sVal =  cos(2*PI*(2*i/wavelength + 2*counter/period));
//      sVal += cos(2*PI*(2*i/wavelength + 1*counter/period));
//      sVal += cos(2*PI*(1*i/wavelength + 0*counter/period));
//      sVal += cos(2*PI*(1*i/wavelength + 1*counter/period));
//     // sVal += cos(2*PI*(3*i/wavelength + 0*counter/period));
//     // sVal += cos(2*PI*(0*i/wavelength + 7*counter/period));
//     // sVal += cos(2*PI*(5*i/wavelength + 9*counter/period));
//     // sVal += cos(2*PI*(6*i/wavelength + 10*counter/period));
//      //sVal += (rand()%100)/70;
//      track2D[counter*trackLength+i]=sVal;
//    }
}


void fft(double *track2D,fftw_complex *outHalf,double *outFull,fftw_plan plan_fwd, int nx, int ny)
{
  int i,j,index,fullIndex,halfIndex;
  int nyh = ny/2+1;
  double modSqr;

//  //print track2D
//  for(i=0;i<nx;i++){
//    printf("row: %2d ",i);
//    for(j=0;j<ny;j++){
//      index = i*trackLength+j;
//      printf("%.2f ",track2D[index]);
//    }
//    printf("\n");
//  }
//  printf("\n\n");


  fftw_execute(plan_fwd);


  //double scaleF = 1;//pow(nx*ny,-1.0);

  for(i=0;i<nx;i++){
    //Build the left ~half of the DFT coef matrix
    for(j=0;j<nyh;j++){
      halfIndex = i*nyh+j;
      fullIndex = i*ny+j;
      modSqr = pow(outHalf[halfIndex][0],2.0)+pow(outHalf[halfIndex][1],2.0);
      outFull[fullIndex] = outFull[fullIndex]+modSqr;
    }

    //Build the right ~half of the DFT coef matrix
    for(j=nyh;j<ny;j++){
      halfIndex = ((nx-i)%nx)*nyh+(ny-j);
      fullIndex = i*ny+j;
      modSqr = pow(outHalf[halfIndex][0],2.0)+pow(outHalf[halfIndex][1],2.0);
      outFull[fullIndex] = outFull[fullIndex]+modSqr;
    }
    //fprintf(outputFile,"\n");
    //printf("\n");
  }
  
    //printf("\n\n");
// PRINT DISTANCE HISTOGRAM
//fprintf(outputFile,"G(r),G_binCount,P(r),P_binCount\n");
//double pofr,gofr; //P(r), G(r)
//for(i=0;i<120;i++){
//  pofr = (float)distNextArray[i]/(counter*numCars);
//  gofr = (float)distAnyArray[i]/(counter*numCars);
//  fprintf(outputFile,"%.8f,%d,%.8f,%d\n",\
//    gofr,distAnyArray[i],pofr,distNextArray[i]);
//}

  //in[0]=0;
  //fftw_execute(plan_bwd);
  //
  //for(i=0;i<nx;i++){
  //  printf("row: %d ",i);
  //  for(j=0;j<ny;j++){
  //    index = i*trackLength+j;
  //    printf("%.2f ",in[index]/(nx*ny));
  //  }
  //  printf("\n");
  //}

}


///////////////////////////////////////////////////////////////
//// Pair Correlation, G(r)
///////////////////////////////////////////////////////////////
void pairCorrelation(int *posArray, int *velArray, int *distArray)
{

  int i,j,dist;

    //Add distances to histogram data
      posArray[numCars] = posArray[0];
      velArray[numCars] = velArray[0];
   
      for (i=0;i<numCars;i++){
        for (j=0;j<numCars;j++){
          dist = posArray[j] - posArray[i];
          while(dist<0){dist = dist + trackLength;}
          distArray[dist]+=1;
        }
      }
}


///////////////////////////////////////////////////////////////
//// Track car diffusion
///////////////////////////////////////////////////////////////
void trackDiffusion(int *posArray, int *velArray)
{
  int i,j,k,step,dist,maxDist;
  int *distArray;

  for (step=0;step<numSteps;step+=sampleRate){

    runSimulation(posArray, velArray);//Will run sampleRate times

    maxDist = 0;
    for (i=0;i<numCars;i++){
      
      dist = posArray[i+1] - posArray[i] - 1;
      while(dist<0){dist = dist + trackLength;}
      if (dist > maxDist)
        maxDist = dist;
    }
    printf("%d,%d\n",step+sampleRate,maxDist);

  } //Finished all steps
}

///////////////////////////////////////////////////////////////
// DISTANCE PROBABITLIY DISTRIBUTION
///////////////////////////////////////////////////////////////
void distProbDistributionHist(int *posArray, int *velArray, int *distArray)
{
  int i,dist;

    //Add distances to histogram data
      posArray[numCars] = posArray[0];
      velArray[numCars] = velArray[0];
   
      for (i=0;i<numCars;i++){
        
        dist = posArray[i+1] - posArray[i];// - 1;
        while(dist<0){dist = dist + trackLength;}
        distArray[dist]+=1;
      }
    

}

///////////////////////////////////////////////////////////////
// PRINT REPRESENTATION OF STATE AS A FUNCTION OF TIME
///////////////////////////////////////////////////////////////
void asciiDisplay(int *posArray, int *velArray)
{
  int i,j,k,step;

  char *track;
  track = malloc((trackLength)*sizeof(char));
  track[trackLength] = '\0';

  for (step=0;step<numSteps;step++){
  
    runSimulation(posArray, velArray);//Will run sampleRate times

    // Print track representation at current state
    for (k=0;k<trackLength;k++)
      track[k] = '.';
    for (k=0;k<numCars;k++)
      track[posArray[k]] = 'X';
    printf("%s\n",track);
    delay(10);
    if(step % 200 == 0) {printf("%d\n",step);}
  }
  free(track);
}


///////////////////////////////////////////////////////////////
// Simulation wrapper
///////////////////////////////////////////////////////////////
void runSimulation(int *posArray, int *velArray)
{
  int i,j;

  for (i=0;i<sampleRate;i++){ //Loop through sampleRate times
    posArray[numCars] = posArray[0];
    velArray[numCars] = velArray[0];
  
    for (j=0;j<numCars;j++){
      updateCar(&posArray[j],&velArray[j],&posArray[j+1]); 
    }
  }
}


///////////////////////////////////////////////////////////////
// Update for one car
///////////////////////////////////////////////////////////////
void updateCar(int *pos1, int *vel1, int *pos2){

int dist,newVel,p;

dist = *pos2 - *pos1 - 1;

while(dist<0){dist = dist + trackLength;}


//Assume speedup if possible
if(*vel1+1 <= maxVel)
  newVel = *vel1+1;
else
  newVel = *vel1;

//Check for available space
if (dist >= newVel){
  *vel1 = newVel; //Set new speed if enough distance
  }
else{ //else set new speed to exactly close gap
  *vel1 = dist;
}

p = 10; //percent chance of slowdown
if (*vel1 > 0){
  if((rand()%100)+1 <= p)
    *vel1 = *vel1 - 1;
}

//update position based on new velocity
*pos1 = (*pos1 + *vel1) % trackLength;

}

float welch(int timeStep)
{
  float n = timeStep;
  float N = fftSamples;

  return 1.0 - pow(((n+1)-.5*(N+1))/(.5*(N+1)),2.0);
}

float Hamm(int timeStep)
{
  float n = timeStep;
  float N = fftSamples;

  return .5*(1-cos(2*PI*n/(N-1)));
}

void delay(int milliseconds)
{
    long pause;
    clock_t now,then;

    pause = milliseconds*(CLOCKS_PER_SEC/1000);
    now = then = clock();
    while( (now-then) < pause )
        now = clock();
}
