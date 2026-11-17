#include "layer.h"
#include "neuron_output.h"
#include "snapshot.h"
#ifndef WIN32
#include <unistd.h>
#endif

#ifdef MACGREGOR_REFERENCE
static const bool referenceModel = true;
#else
static const bool referenceModel = false;
#endif

void Layer::allocateConnections() {
    ConnectionMatrix = new std::vector<Synapse>*[rows]();
    for (int i = 0; i < rows; ++i)
        ConnectionMatrix[i] = new std::vector<Synapse>[clms]();
}

Layer::~Layer() {
    for (int i = 0; i < rows; ++i) {
        if (ConnectionMatrix) delete[] ConnectionMatrix[i];
        if (NeuronLayer) delete[] NeuronLayer[i];
    }
    delete[] ConnectionMatrix;
    delete[] NeuronLayer;
}

int Layer::simulate(Param& param, const RunOptions& options, Topology topology) try {
#ifdef WIN32
    if (options.snapshotSteps || !options.resume.empty())
        throw std::runtime_error("Restart snapshots require the POSIX drand48 generator");
#endif
    int firstStep = 0;
    EpochState epoch(options.epochSteps, firstStep, rows * clms);
    std::string prefix;
    if (options.resume.empty()) {
        addNeurons(&param);
        addSynapse(0, &param);
    } else {
        NetworkSnapshot snapshot(options.resume, true);
        snapshot.header(topology, referenceModel, rows, clms, firstStep, param);
        if (firstStep >= options.steps)
            throw std::runtime_error("--steps must exceed the snapshot's completed step count");
        addNeurons(&param);
        allocateConnections();
        snapshot.state(rows, clms, NeuronLayer, ConnectionMatrix, firstStep, epoch);
        if (epoch.interval != options.epochSteps)
            epoch = EpochState(options.epochSteps, firstStep, rows * clms);
        prefix = "from-" + std::to_string(firstStep) + "-";
    }

    writeEdges(prefix, rows, clms, ConnectionMatrix);
    RunOutput output(prefix, rows, clms, 2.0, epoch);
    NeuronTraces traces(perNeuronOutput, rows * clms, prefix);
    const int traceStart = topology == Topology::smallworld ? 0 : 5001;
    const int traceEnd = topology == Topology::scalefree ? 10000 :
                         topology == Topology::uniform ? 65000 : INT_MAX;
    for (int step = firstStep; step < options.steps; ++step) {
        // Preserve the original draw for a potential stimulus (the applied current is zero).
        (void) std::rand();
        if (perNeuronOutput && step >= traceStart && step < traceEnd)
            traces.writeStep(step, rows, clms, NeuronLayer, ConnectionMatrix);
        output.recordSpikes(step, NeuronLayer);

#pragma omp parallel for
        for (int row = 0; row < rows; ++row)
            for (int column = 0; column < clms; ++column)
                NeuronLayer[row][column].glutamaterelease(ConnectionMatrix[row][column], NeuronLayer, step);
#pragma omp parallel for
        for (int row = 0; row < rows; ++row)
            for (int column = 0; column < clms; ++column)
                NeuronLayer[row][column].computePotential(ConnectionMatrix[row][column], 0, step);

        output.recordField(step, NeuronLayer, referenceModel);
        int completed = step + 1;
        if (epoch.interval && completed - epoch.firstStep == epoch.interval)
            output.finishEpoch(completed);
        if (options.snapshotSteps && completed % options.snapshotSteps == 0) {
            const std::string path = prefix + "snapshot-" + std::to_string(completed) + ".bin";
            // Publish a checkpoint only after its corresponding observations reach disk buffers.
            if (!std::cout.flush() || std::fflush(nullptr) != 0)
                throw std::runtime_error("Cannot flush recordings before checkpoint");
            std::string temporary = path + ".tmp";
#ifndef WIN32
            temporary += "." + std::to_string(getpid());
#endif
            try {
                NetworkSnapshot snapshot(temporary, false);
                snapshot.header(topology, referenceModel, rows, clms, completed, param);
                snapshot.state(rows, clms, NeuronLayer, ConnectionMatrix, completed, epoch);
                if (std::rename(temporary.c_str(), path.c_str()) != 0)
                    throw std::runtime_error(path + ": cannot publish snapshot");
            } catch (...) {
                std::remove(temporary.c_str());
                throw;
            }
        }
    }
    if (epoch.interval && epoch.firstStep < options.steps)
        output.finishEpoch(options.steps);
    return EXIT_SUCCESS;
} catch (const std::exception& error) {
    std::fprintf(stderr, "%s\n", error.what());
    return EXIT_FAILURE;
}
