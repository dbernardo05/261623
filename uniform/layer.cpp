
#include "layer.h"
#include "../neuron_output.h"

double distribution_max;

Layer::Layer(int row, int clm, bool perNeuronOutput)
    : rows(row), clms(clm), perNeuronOutput(perNeuronOutput) {}
//add Neurons to grid
int Layer::addNeurons(Param* param){
	int i, j;
	double mean, stddev;
	distribution_max = 0;
	mean = 50; stddev = 10.0;
	default_random_engine generator;
	normal_distribution<double> distribution(mean,stddev);
	//The initial plan was to pass the param pointer to the constructor. However, that would have involved initializing the neurons in a loop, as C++ does not have a method to initialize an array created with new without passing param to each instance in a loop.
	NeuronLayer = new Neuron*[rows]();
	for (i = 0; i < rows; i++){ 
		NeuronLayer[i] = new Neuron[clms]();
	}
	for(i = 0; i < rows; i++){
		for(j =0; j < clms; j++){
			NeuronLayer[i][j].row = i;
			NeuronLayer[i][j].clm = j;		
			if((i%5 == 0) && (j%5 == 0)){
				NeuronLayer[i][j].excitatory = false;
				NeuronLayer[i][j].Inhibitory();
				NeuronLayer[i][j].Tgk = param->InTgk;
				NeuronLayer[i][j].Tth = param->InTth;
				NeuronLayer[i][j].B   = param->InB;
				NeuronLayer[i][j].Te1 = param->InTe1;
				NeuronLayer[i][j].Te2 = param->InTe2;
				NeuronLayer[i][j].Ti1 = param->InTi1;
				NeuronLayer[i][j].Ti2 = param->InTi2;
				//cout << i << " " << j << endl;
			}
			else{
				NeuronLayer[i][j].excitatory = true;
				NeuronLayer[i][j].Tgk = param->ExTgk;
				NeuronLayer[i][j].Tth = param->ExTth;
				NeuronLayer[i][j].B   = param->ExB;
				NeuronLayer[i][j].Te1 = param->ExTe1;
				NeuronLayer[i][j].Te2 = param->ExTe2;
				NeuronLayer[i][j].Ti1 = param->ExTi1;
				NeuronLayer[i][j].Ti2 = param->ExTi2;
				NeuronLayer[i][j].TauRelease	= param->TauRelease;
				NeuronLayer[i][j].TauReplenish	= param->TauReplenish;
				NeuronLayer[i][j].SpontRelease	= param->SpontRelease;
				NeuronLayer[i][j].glutamate_pr = distribution(generator);
				if(distribution_max < NeuronLayer[i][j].glutamate_pr)
					distribution_max = NeuronLayer[i][j].glutamate_pr;
			}
		}
	}
	return 0;
}
//add synaptic connections
int Layer::addSynapse(int no_wts, Param* param){
	int i, j, k, m, n, count, connect;
	double d = param->Neighborhood, value, value_uniform,p, in_Glut, wt;
	const auto glutamateMap = openOutputFile("glutamatemap.txt", perNeuronOutput);
	const auto initialGlutamateMap = openOutputFile("initialglutamatemap.txt", perNeuronOutput);
	const auto weightMap = openOutputFile("weightmap.txt", perNeuronOutput);
	vector<double> weightRow(perNeuronOutput ? rows * clms : 0);
	allocateConnections();
	
//#pragma omp parallel for private(i,j,k,p,value,count,m,n)
	for(i = 0; i < rows; i++){
		for(j = 0; j < clms; j++){
			count = 0;
			in_Glut = 0;
			for(m = 0; m < rows; m++){
				for(n = 0; n < clms; n++){
					//initial connection criteria: distance between 2 neurons
					//no self feedback
					//(sqrt((i - m)*(i - m) + (j - n)*(j - n)) < d) &&
					value = drand48();
					value_uniform = drand48();
					//(value_uniform < 0.625) && 
					if((value < exp(-((i - m)*(i - m))/5 -((j - n)*(j - n))/5)) &&  (!(i == m && j == n)) && (!(!NeuronLayer[i][j].excitatory && !NeuronLayer[m][n].excitatory)) ){
					ConnectionMatrix[i][j].push_back(Synapse());
					ConnectionMatrix[i][j][count].row = m;
					ConnectionMatrix[i][j][count].clm = n;
					value = drand48();
					ConnectionMatrix[i][j][count].weight = value;
					if(NeuronLayer[m][n].excitatory){
						ConnectionMatrix[i][j][count].excitatory = true;
						ConnectionMatrix[i][j][count].NRmax = (param->NRmax)*value;
						ConnectionMatrix[i][j][count].tau	= 0;
					}
					else{
						ConnectionMatrix[i][j][count].excitatory = false;
						ConnectionMatrix[i][j][count].NRmax = 0;
						ConnectionMatrix[i][j][count].tau	= 0;
						ConnectionMatrix[i][j][count].NR    = 0;
					}
					if (perNeuronOutput) {
						weightRow[m * clms + n] = value;
						in_Glut += ConnectionMatrix[i][j][count].NRmax;
					}
					count = count+1;
					//fprintf(cFile,"%f ",value);
					}					
				}
			}		
  			//cout << count << " " ;
			if (perNeuronOutput) {
				fprintf(glutamateMap.get(), "%f ", in_Glut);
				writeWeightRow(weightMap.get(), weightRow);
			}
		}
		if (perNeuronOutput)
			fputc('\n', glutamateMap.get());
	}
//exit(0);	
//#pragma omp parallel for private(i,j)
	for(i = 0; i < rows; i++){
		for(j = 0; j < clms; j++){
			NeuronLayer[i][j].glutamateinitialize(ConnectionMatrix[i][j], NeuronLayer);
		}
	}
	//exit(0);	
	for(i = 0; i < rows; i++){
		for(j = 0; j < clms; j++){
			in_Glut = 0;
			for(count = 0; count < ConnectionMatrix[i][j].size(); count++){
				if (perNeuronOutput)
					in_Glut += ConnectionMatrix[i][j][count].NR;
				if(ConnectionMatrix[i][j][count].excitatory)
					NeuronLayer[ConnectionMatrix[i][j][count].row][ConnectionMatrix[i][j][count].clm].glutamate += ConnectionMatrix[i][j][count].NR;
			}
			if (perNeuronOutput)
				fprintf(initialGlutamateMap.get(), "%f ", in_Glut);
		}
		if (perNeuronOutput)
			fputc('\n', initialGlutamateMap.get());
	}
	//exit(0);
	return 0;
}
