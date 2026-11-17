#pragma once

#include "output_options.h"
#include <cerrno>
#include <climits>
#include <string>

struct RunOptions {
    bool perNeuronOutput = false;
    int steps = 60000;
    int epochSteps = 0;
    int snapshotSteps = 0;
    std::string resume;
};

enum class Topology { uniform = 1, smallworld, scalefree };

inline RunOptions parseRunOptions(int argc, char* argv[]) {
    const char* usage =
        "Usage: %s [--per-neuron-output] [--steps N] [--epoch-steps N]\n"
        "          [--snapshot-steps N] [--resume FILE]\n"
        "  --per-neuron-output  Write the dense maps and traces (off by default).\n"
        "  --steps N            Stop after step N-1 (default 60000).\n"
        "  --epoch-steps N      Summarize each N steps, including the final partial epoch.\n"
        "  --snapshot-steps N   Save full network state every N completed steps.\n"
        "  --resume FILE        Continue a snapshot with its original step numbering.\n";
    RunOptions options;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--help" && argc == 2) {
            std::printf(usage, argv[0]);
            std::exit(EXIT_SUCCESS);
        }
        if (arg == "--per-neuron-output" && !options.perNeuronOutput) {
            options.perNeuronOutput = true;
            continue;
        }
        if (arg == "--resume" && i + 1 < argc && argv[i+1][0] && argv[i+1][0] != '-') {
            options.resume = argv[++i];
            continue;
        }
        if ((arg == "--steps" || arg == "--epoch-steps" || arg == "--snapshot-steps") && i + 1 < argc) {
            const char* word = argv[++i];
            char* end;
            errno = 0;
            const long n = std::strtol(word, &end, 10);
            if (!errno && *word && !*end && n > 0 && n <= INT_MAX) {
                if (arg == "--steps") options.steps = static_cast<int>(n);
                if (arg == "--epoch-steps") options.epochSteps = static_cast<int>(n);
                if (arg == "--snapshot-steps") options.snapshotSteps = static_cast<int>(n);
                continue;
            }
        }
        std::fprintf(stderr, "Invalid option or value for %s\n", arg.c_str());
        std::fprintf(stderr, usage, argv[0]);
        std::exit(EXIT_FAILURE);
    }
    return options;
}
