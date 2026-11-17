#pragma once

#include "output_options.h"
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

struct FieldPotential {
    double legacy = 0;
    double fixed = 0;
};

class Electrodes {
    int rows, columns;
    std::vector<FieldPotential> weights;

public:
    Electrodes(int rows, int columns, double sigma)
        : rows(rows), columns(columns), weights(rows * columns) {
        const double scale = 1.0 / std::sqrt(2 * std::acos(-1.0) * sigma * sigma);
        for (int row = 0; row < rows; ++row) {
            for (int column = 0; column < columns; ++column) {
                const int dy = rows / 2 - row, dx = columns / 2 - column;
                const double r2 = dy * dy + dx * dx;
                auto& w = weights[row * columns + column];
                // The legacy expression multiplies by sigma twice after dividing by two.
                w.legacy = scale * std::exp(-r2 / 2.0 * sigma * sigma);
                w.fixed = scale * std::exp(-r2 / (2.0 * sigma * sigma));
            }
        }
    }

    template<class Neuron>
    FieldPotential sample(const Neuron* const* neurons) const {
        FieldPotential field;
        for (int row = 0; row < rows; ++row) {
            for (int column = 0; column < columns; ++column) {
                const auto& w = weights[row * columns + column];
                const double e = neurons[row][column].E;
                field.legacy += w.legacy * e;
                field.fixed += w.fixed * e;
            }
        }
        return field;
    }
};

template<class Synapse>
void writeEdges(const std::string& prefix, int rows, int columns,
                const std::vector<Synapse>* const* connections) {
    const auto file = openOutputFile(prefix + "edges.tsv");
    std::fputs("pre_neuron\tpost_neuron\tweight\texcitatory\n", file.get());
    for (int row = 0; row < rows; ++row) {
        for (int column = 0; column < columns; ++column) {
            for (const auto& edge : connections[row][column])
                std::fprintf(file.get(), "%d\t%d\t%.17g\t%d\n", edge.row * columns + edge.clm,
                             row * columns + column, edge.weight, int(edge.excitatory));
        }
    }
}

struct EpochState {
    int interval, firstStep;
    std::int64_t excitatory = 0, inhibitory = 0;
    FieldPotential fields;
    std::vector<int> spikes;
    std::vector<double> glutamate;

    EpochState(int interval, int firstStep, int cells)
        : interval(interval), firstStep(firstStep), spikes(interval ? cells : 0), glutamate(spikes.size()) {}
};

class RunOutput {
    int rows, columns;
    OutputFile counts, events, fields, epochs, neuronEpochs;
    Electrodes electrodes;
    EpochState& epoch;

public:
    RunOutput(const std::string& prefix, int rows, int columns, double sigma, EpochState& epoch)
        : rows(rows), columns(columns), counts(openOutputFile(prefix + "spikes.txt")),
          events(openOutputFile(prefix + "spike_events.tsv")), fields(openOutputFile(prefix + "field_potentials.tsv")),
          epochs(openOutputFile(prefix + "epochs.tsv", epoch.interval != 0)),
          neuronEpochs(openOutputFile(prefix + "neuron_epochs.tsv", epoch.interval != 0)),
          electrodes(rows, columns, sigma), epoch(epoch) {
        std::fputs("step\tneuron\n", events.get());
        std::fputs("step\tlegacy\tfixed\n", fields.get());
        if (epochs) {
            std::fputs("first_step\tend_step\texcitatory_spikes\tinhibitory_spikes\tmean_legacy\tmean_fixed\n", epochs.get());
            std::fputs("first_step\tend_step\tneuron\tspikes\tmean_received_glutamate\n", neuronEpochs.get());
        }
    }

    // Match spikes.txt and the dense spike patterns, which sample S before the update.
    template<class Neuron>
    void recordSpikes(int step, const Neuron* const* neurons) {
        int excitatory = 0, inhibitory = 0;
        for (int row = 0; row < rows; ++row) {
            for (int column = 0; column < columns; ++column) {
                const auto& neuron = neurons[row][column];
                if (!neuron.S) continue;
                const int id = row * columns + column;
                std::fprintf(events.get(), "%d\t%d\n", step, id);
                if (neuron.excitatory) ++excitatory;
                else ++inhibitory;
                if (epochs) ++epoch.spikes[id];
            }
        }
        std::fprintf(counts.get(), "%d %d %d\n", step, excitatory, inhibitory);
        if (epochs) {
            epoch.excitatory += excitatory;
            epoch.inhibitory += inhibitory;
        }
    }

    template<class Neuron>
    void recordField(int step, const Neuron* const* neurons, bool legacyStdout) {
        const auto field = electrodes.sample(neurons);
        std::fprintf(fields.get(), "%d\t%.17g\t%.17g\n", step, field.legacy, field.fixed);
        std::cout << (legacyStdout ? field.legacy : field.fixed) << ' ' << step << '\n';
        if (epochs) {
            epoch.fields.legacy += field.legacy;
            epoch.fields.fixed += field.fixed;
            for (int row = 0; row < rows; ++row)
                for (int column = 0; column < columns; ++column)
                    epoch.glutamate[row * columns + column] += neurons[row][column].glutamate;
        }
    }

    void finishEpoch(int endStep) {
        const int length = endStep - epoch.firstStep;
        std::fprintf(epochs.get(), "%d\t%d\t%lld\t%lld\t%.17g\t%.17g\n", epoch.firstStep, endStep,
                     static_cast<long long>(epoch.excitatory), static_cast<long long>(epoch.inhibitory),
                     epoch.fields.legacy / length, epoch.fields.fixed / length);
        for (size_t id = 0; id < epoch.spikes.size(); ++id) {
            std::fprintf(neuronEpochs.get(), "%d\t%d\t%zu\t%d\t%.17g\n", epoch.firstStep, endStep, id,
                         epoch.spikes[id], epoch.glutamate[id] / length);
            epoch.spikes[id] = 0;
            epoch.glutamate[id] = 0;
        }
        epoch.firstStep = endStep;
        epoch.excitatory = epoch.inhibitory = 0;
        epoch.fields = FieldPotential();
    }
};
