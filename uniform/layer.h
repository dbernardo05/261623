#include "neuron.h"
#include "../run_options.h"


class Layer{
	int rows, clms;
	bool perNeuronOutput;
	//NeuronLayer is a matrix of neurons
	Neuron **NeuronLayer = nullptr;
	//ConnectionMatrix is a matrix of vectors; each vector contains the incoming synapses to the neuron at that location
	vector<Synapse> **ConnectionMatrix = nullptr;
	void allocateConnections();
public:
	Layer(int row, int clm, bool perNeuronOutput = false);
	int addNeurons(Param* param);
	int addSynapse(int no_wts, Param* param);
	~Layer();
	Layer(const Layer&) = delete;
	Layer& operator=(const Layer&) = delete;
	int simulate(Param& param, const RunOptions& options, Topology topology);
};
