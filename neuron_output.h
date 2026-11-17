#pragma once

#include "output_options.h"
#include <vector>

// Incoming weights for one postsynaptic neuron, in presynaptic grid order.
inline void writeWeightRow(std::FILE* file, std::vector<double>& weights) {
    for (double& weight : weights) {
        std::fprintf(file, "%f ", weight);
        weight = 0;
    }
    std::fputc('\n', file);
}

class NeuronTraces {
    OutputFile potential, releasedGlutamate, releasableGlutamate, availableGlutamate;
    OutputFile gaba, spikePatterns, kConductance;
    std::vector<double> availableByNeuron;

public:
    NeuronTraces(bool enabled, int neuronCount, const std::string& prefix)
        : potential(openOutputFile(prefix + "potential.txt", enabled)),
          releasedGlutamate(openOutputFile(prefix + "releasedglutamate.txt", enabled)),
          releasableGlutamate(openOutputFile(prefix + "releasableglutamate.txt", enabled)),
          availableGlutamate(openOutputFile(prefix + "availableglutamate.txt", enabled)),
          gaba(openOutputFile(prefix + "gaba.txt", enabled)),
          spikePatterns(openOutputFile(prefix + "spikepatterns.txt", enabled)),
          kConductance(openOutputFile(prefix + "kconductance.txt", enabled)),
          availableByNeuron(enabled ? neuronCount : 0) {}

    template<class Neuron, class Synapse>
    void writeStep(int step, int rows, int columns, const Neuron* const* neurons,
                   const std::vector<Synapse>* const* connections) {
        for (int row = 0; row < rows; ++row) {
            for (int column = 0; column < columns; ++column) {
                double releasable = 0;
                for (const auto& synapse : connections[row][column]) {
                    releasable += synapse.NR;
                    if (synapse.excitatory)
                        availableByNeuron[synapse.row * columns + synapse.clm] += synapse.NR;
                }
                const auto& neuron = neurons[row][column];
                std::fprintf(potential.get(), "%f ", neuron.E);
                std::fprintf(releasedGlutamate.get(), "%f ", neuron.glutamate);
                std::fprintf(releasableGlutamate.get(), "%f ", releasable);
                // The upstream neuron first sets GABA during the update at step zero.
                std::fprintf(gaba.get(), "%f ", step == 0 ? 0.0 : neuron.gaba);
                std::fprintf(spikePatterns.get(), "%d ", neuron.S);
                std::fprintf(kConductance.get(), "%f ", neuron.Gk);
            }
        }
        for (double& available : availableByNeuron) {
            std::fprintf(availableGlutamate.get(), "%f ", available);
            available = 0;
        }
        for (std::FILE* file : {potential.get(), releasedGlutamate.get(), releasableGlutamate.get(),
                                availableGlutamate.get(), gaba.get(), spikePatterns.get(), kConductance.get()})
            std::fputc('\n', file);
    }
};
