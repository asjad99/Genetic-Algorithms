#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <cstdio>

using namespace std;


const int cPopSize = 100; // must be an even number
const int cNumGens = 150;
const double cMutationRate = 0.01;
const double cCrossoverRate = 0.70;
int TOTAL_CITIES = 0;

float **lookup_table;


void InitPopulation(int ***CrntPop,int ***NextPop,float **Fitness,int **BestMenber);

void CreateLookUpTable(int *xcord_array, int *ycord_array);
float calc_distance(int dX0, int dY0, int dX1, int dY1);
int roulette_wheel_select(float *Fitness);
void Crossover(int *P1,int *P2,int *C1,int *C2);
void Copy(int *P1,int *P2,int *C1,int *C2);

void FreeMem(int **CrntPop,int **NextPop,float *Fitness,int *BestMember);
float **Aloc2DAry(int m,int n);
void Free2DAry(float **Ary2D,int n);

void Mutate(int ***CrntPop);

long random_gen(long max);
double Rand01();
int RandInt(int n);

int main(int argc,char *argv[]){

	/*------------------------------
	GA approach:

		1) Initialize Population
		2) Evaluate the fitness of the population
		3) Generate new population by Evolution(Crossover and Mutation)  
		
	//------------------------------*/


	//read a data file containing the locations of the cities into a dynamically allocated array
	//-------------------------

	ifstream fin;
	char FName[20];
	int total_cities = 0;
	int Tmp_input = 0;
	
	cout<<"Enter data filename: ";
    cin>>FName; cin.ignore();
  
  	//find out the number of cities and allocate the array. 
  	fin.open(FName);
  	
  	if(!fin.good()){cout<<"File not found!\n";exit(1);}

  	fin>>Tmp_input;     //Number of outputs 

  	TOTAL_CITIES = Tmp_input;
  	cout << "Total Cities" << TOTAL_CITIES <<'\n';

  	int *xcord_array = new int[TOTAL_CITIES];
	int *ycord_array =  new int[TOTAL_CITIES];

  	for(int index=0;index<TOTAL_CITIES;index++){
		fin>>xcord_array[index];
		fin>>ycord_array[index];
	}
  	//for testing...TOTAL_CITIES = 10;

  	fin.close();

  	//for(int index=0;index<TOTAL_CITIES;index++)
  	//	cout << xcord_array[index] << "    "<<ycord_array[index]  <<"\n";
  	

  	//---------------------(Initialization)----------------
  	
  	
  	int **CrntPop, **NextPop; // the crnt & next population lives here
  	float *Fitness;
  	int BestFitness=0, *BestMember; // fitness vars
  	int i,j = 0;

  	InitPopulation (&CrntPop,&NextPop,&Fitness,&BestMember);

  	CreateLookUpTable(xcord_array,ycord_array);

    //---------------------(Evaluate Fitness)----------------
    //each iteration of our main loop represents a new generation. 

    for(int Gen=0;Gen<cNumGens;Gen++){	

    	
  		int city_number1 = 0;
  		int city_number2 = 0;
  		double fitness_temp = 0; 

  		//calculate the route distances(fitness of each member in GA terms) using the lookup table
  		for(i=0;i<cPopSize;i++){
  			//cout <<"\n---------------\n";
    		
    		for(j=0;j<TOTAL_CITIES;j++){

    			if (j == TOTAL_CITIES-1){
    				city_number1 = CrntPop[i][j];
    				city_number2 = CrntPop[i][0];
    			}
    			else{
    				city_number1 = CrntPop[i][j];
    				city_number2 = CrntPop[i][j+1];
    			}
    			
    			if (city_number2 >TOTAL_CITIES){
    				city_number2 = RandInt(TOTAL_CITIES);
    			}


    			if (city_number1 >TOTAL_CITIES){
    				city_number1 = RandInt(TOTAL_CITIES);
    			}

    			//cout <<"Gen:" << Gen <<"pop member" <<i << "index:" << j << ":" << "city_number1:" << city_number1 << "city_number2:" << city_number2  <<"\n";
      
      			fitness_temp += lookup_table[city_number1][city_number2];
      			
    		}	
    		Fitness[i] =  1/fitness_temp * 1000000;

    		if(BestFitness < Fitness[i]){ // save best member
       			 BestFitness=Fitness[i];
       			  for(int j=0;j<TOTAL_CITIES;j++)BestMember[j]=CrntPop[i][j];
       		}

    		//cout << "fitness:" << Fitness[i] <<'\n'; 
  		}
    
  	//---------------------(Produce the next population)----------------

  		for(i=0;i<cPopSize;i+=2){

  		 	int Parent1=roulette_wheel_select(Fitness);
      		int Parent2=roulette_wheel_select(Fitness);

      		//cout << "-----------------\n";
      		//cout << Parent1 << Parent2 << "\n";

      		//TODO: We then check if crossover occurs (it should occur 70 percent of the time)
      		if(cCrossoverRate>Rand01()){
      			Crossover(CrntPop[Parent1],CrntPop[Parent2],NextPop[i],NextPop[i+1]);
	  		}
     	    else
       			Copy(CrntPop[Parent1],CrntPop[Parent2],NextPop[i],NextPop[i+1]);
     	

     		if(cMutationRate<Rand01())Mutate(&NextPop);
  		 }
  		 	int **Tmp  = CrntPop; 
  			CrntPop = NextPop; 
    		NextPop = Tmp;

  			cout<<setw(3)<<Gen<<':'<<setw(5)<<BestFitness<<endl;
  	}

  	cout<<"Best Fitness:" << BestFitness <<"!\n";
  	
    cout<<"Best Individual: ";
  	for(i=0;i<TOTAL_CITIES;i++)
  		cout<<BestMember[i] <<'\n';cout<<endl;
  
  	//Free all used memory
  	 FreeMem(CrntPop,NextPop,Fitness,BestMember);

  	 //free 2d arrays,lookup tables etc.

  	 delete [] xcord_array;
  	 delete [] ycord_array;
  	 //delete [] route_distance;

  	 Free2DAry(lookup_table,TOTAL_CITIES);

} //end of main

//----------------------------------------------------------------------
void  InitPopulation (int ***CrntPop,int ***NextPop,float **Fitness,int **BestMember){

  //1. declare two arrays of population size(100) length. one is for current population, one is for next population
  //2. each item holds the address to the  

  int i, j;
  //srand(Seed);
  *CrntPop = new int*[cPopSize];
  *NextPop = new int*[cPopSize];
  
  for(i=0;i<cPopSize;i++){
    (*CrntPop)[i] = new int[TOTAL_CITIES]; 
    (*NextPop)[i] = new int[TOTAL_CITIES];
  }
  *Fitness    = new float[cPopSize]; //fitness var points to an int array of size 100. 
  *BestMember = new int[TOTAL_CITIES];  //Best member var points to an int array of size 80.
  
  if(Fitness==NULL||BestMember==NULL)exit(1);
 
  for(i=0;i<cPopSize;i++){
    for(j=0;j<TOTAL_CITIES;j++){
      (*CrntPop)[i][j] = j;
    }
  }

  for(i=0;i<cPopSize;i++){
    for(j=0;j<TOTAL_CITIES;j++){
      int temp_swap = 0;
      int random_num = RandInt(TOTAL_CITIES);

     
      //cout << (*CrntPop)[i][j];
      temp_swap = (*CrntPop)[i][j];

      (*CrntPop)[i][j] = (*CrntPop)[i][random_num];
      (*CrntPop)[i][random_num] = temp_swap;

    }
  }

  for(i=0;i<cPopSize;i++){
    for(j=0;j<TOTAL_CITIES;j++){
      int temp_swap = 0;
      int random_num = RandInt(TOTAL_CITIES);
     
      temp_swap = (*NextPop)[i][j];

      (*NextPop)[i][j] = (*NextPop)[i][random_num];
      (*NextPop)[i][random_num] = temp_swap;

    }
  }
    	
 }
//----------------------------------------------------------------------
//n x n lookup table
void CreateLookUpTable(int *xcord_array, int *ycord_array){

	//for(int index=0;index<TOTAL_CITIES;index++)
  	//	cout << xcord_array[index] << "    "<<ycord_array[index]  <<"\n";

  	//--------create 2d lookup table------

  	lookup_table   = Aloc2DAry(TOTAL_CITIES,TOTAL_CITIES);
  	
  	//--------populate the array------
  	//calcuate the distance between city 1 to city 1,2,3,4,5, city 2 to 1,2,3,4,5,6 and so on (1 to 1 and 2 to 2 would be set to zero)

  	for(int i=0;i<TOTAL_CITIES;i++){

   	 	for(int j=0;j<TOTAL_CITIES;j++){
   	 		lookup_table[i][j]	= calc_distance(xcord_array[i],ycord_array[i],xcord_array[j],ycord_array[j]) ; 
	
   	 	}
  	}    
}
//calculate distance between two cities
//----------------------------------------------------
float calc_distance(int dX0, int dY0, int dX1, int dY1) {
    return sqrt((dX1 - dX0)*(dX1 - dX0) + (dY1 - dY0)*(dY1 - dY0));
}

//----------------------------------------------------
int roulette_wheel_select(float *Fitness){
	//Copy fitnesses into a temp array
		//Subtract the minimum pop fitness from each fitness to normalise them
		//Subtract each normalised fitness from the normalized maximum pop fitness to get the scaled fitnesses Generate a random number between 0 and the sum of all the scaled fitnesses

		int i;
		int j;
		float *relative_fitness = new float[TOTAL_CITIES];
		float *prob = new float[TOTAL_CITIES];
		
		double total_fitness = 0;
		
		//calc. total fitness
		for (i=0;i<cPopSize;i++){
			total_fitness += Fitness[i];
			//if(TempSum>RandInt(cPopSize)) return i; 
		}

		//calc. relative fitness of each member
		for (i=0;i<cPopSize;i++){
			relative_fitness[i] = Fitness[i]/total_fitness;
			//cout << Fitness[i] <<'\n';
			//cout <<"relative:"<< relative_fitness[i] <<'\n';
		}

		// Generate probability intervals for each individual
    	for (i=0;i<cPopSize;i++){

			for(j=0;j<=i;j++) { 
				prob[i] += relative_fitness[j];
			}
			//cout <<"fitness:" << Fitness[i] <<"prob:"<< prob[i] <<'\n';
		}


		double r = ((double) rand() / (RAND_MAX));

		for (i=0;i<cPopSize;i++){
			if (r <= prob[i]){
				return i;
			}
		}    

}

void Crossover(int *P1,int *P2,int *C1,int *C2){

	//uses moon crossover technique to generate offspring 
	int i=0;
	int j=0;
	int rand_point = RandInt(TOTAL_CITIES);

	//pick a random point in parent 1 array (CrntPop[Parent1])
	//copy all points to the left of parent 1 array to the child 1 array

	for (i=0; i<rand_point;i++){
		C1[i]=P1[i];
	}

	//start copying elements from parent2, leaving ones that already exist(to avoid dublicates)
	int parent2_index = 0;
	int counter_check = 0;
	for (i=rand_point; i<TOTAL_CITIES;i++) {

		for (j=0; j<rand_point;j++){
			
		 	//dublicate found..pick next element from parent 2
		 	if (C1[j] == P2[parent2_index]){
		 			
		 			if (parent2_index > 99){
						parent2_index = 0;
					}
					else{
						parent2_index++;
					}
					//start over 
		 			j=0;
		 			
		 	}
		 	//cout << "j:  " <<j << "   parent_index 2:"<< parent2_index <<"\n";
		}
		//cout << "copying...";
		
		if (parent2_index < 100){
			C1[i]=P2[parent2_index];
			parent2_index++;
		}
		else{
			parent2_index = 0;
			
		}

	}

	int rand_point_2 = RandInt(TOTAL_CITIES);

	//pick a random point in parent 1 array (CrntPop[Parent1])
	//copy all points to the left of parent 1 array to the child 1 array

	for (i=0; i<rand_point_2;i++){
		C2[i]=P2[i];
	}

	//start copying elements from parent2, leaving ones that already exist(to avoid dublicates)
	int parent1_index = 0;
	for (i=rand_point_2; i<TOTAL_CITIES;i++) {
		for (j=0; j<rand_point_2;j++){
		 	//dublicate found..pick next one
		 	if (C1[j]==P2[parent1_index]){
		 		
		 			if (parent1_index > 99){
						parent1_index = 0;
					}
					else{
						parent1_index++;
					}
		 			j=0;
		 	}
		 	//cout << "j:  " <<j << "   parent_index 1:"<< parent1_index <<"\n";
		
		}
		if (parent1_index < 100){
			C2[i]=P1[parent1_index];
			parent1_index++;
		}
		else{
			parent1_index = 0;
			
		}
	}
}

//----------------------------------------------------
int RandInt(int n){ // 0..n-1
  return int( rand()/(double(RAND_MAX)+1) * n );
}


void Copy(int *P1,int *P2,int *C1,int *C2){
  for(int i=0;i<TOTAL_CITIES;i++){
    C1[i]=P1[i]; C2[i]=P2[i];
  }
}

void Mutate(int ***CrntPop) {
int i =0;
int j =0;
  for(i=0;i<cPopSize;i++){
    for(j=0;j<TOTAL_CITIES;j++){
      int temp_swap = 0;
      int random_num = RandInt(TOTAL_CITIES);

      temp_swap = (*CrntPop)[i][j];

      (*CrntPop)[i][j] = (*CrntPop)[i][random_num];
      (*CrntPop)[i][random_num] = temp_swap;

    }
}

}

//---------------------------------------------------------------------


float **Aloc2DAry(int m,int n){
//Allocates memory for 2D array
  float **Ary2D = new float*[m];
  if(Ary2D==NULL){cout<<"No memory!\n";exit(1);}
  for(int i=0;i<m;i++){
	 Ary2D[i] = new float[n];
	 if(Ary2D[i]==NULL){cout<<"No memory!\n";exit(1);}
  }
  return Ary2D;
}

void Free2DAry(float **Ary2D,int n){
//Frees memory in 2D array
  for(int i=0;i<n;i++)
	 delete [] Ary2D[i];
  delete [] Ary2D;
}


//----------------------------------------------------------------------
void FreeMem(int **CrntPop,int **NextPop,float *Fitness,int *BestMenber){
  for(int i=0;i<cPopSize;i++){
    delete[]CrntPop[i];
    delete[]NextPop[i];
  }
  delete CrntPop;
  delete NextPop;
  delete Fitness;
  delete BestMenber;
}

double Rand01(){ // 0..1
  return(rand()/(double)(RAND_MAX));
}


//----------------------------------------------------------------------

// Assumes 0 <= max <= RAND_MAX
// Returns in the half-open interval [0, max]
long random_gen(long max) {
  unsigned long
    // max <= RAND_MAX < ULONG_MAX, so this is okay.
    num_bins = (unsigned long) max + 1,
    num_rand = (unsigned long) RAND_MAX + 1,
    bin_size = num_rand / num_bins,
    defect   = num_rand % num_bins;

  long x;
  do {
   x = random();
  }
  // This is carefully written not to overflow
  while (num_rand - defect <= (unsigned long)x);

  // Truncated division is intentional
  return x/bin_size;
}