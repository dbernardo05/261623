#pragma once

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>

inline bool perNeuronOutputRequested(int argc, char* argv[]) {
    if (argc == 1)
        return false;
    if (argc == 2 && std::strcmp(argv[1], "--per-neuron-output") == 0)
        return true;

    const bool help = argc == 2 && std::strcmp(argv[1], "--help") == 0;
    std::fprintf(help ? stdout : stderr,
                 "Usage: %s [--per-neuron-output]\n"
                 "  --per-neuron-output  Write individual neuron data (off by default).\n",
                 argv[0]);
    std::exit(help ? EXIT_SUCCESS : EXIT_FAILURE);
}

struct CloseOutputFile {
    std::string path;
    void operator()(std::FILE* file) const {
        const bool failed = std::ferror(file) != 0;
        if (std::fclose(file) != 0 || failed) {
            std::fprintf(stderr, "Could not finish writing %s\n", path.c_str());
            std::exit(EXIT_FAILURE);
        }
    }
};
using OutputFile = std::unique_ptr<std::FILE, CloseOutputFile>;

inline OutputFile openOutputFile(const std::string& path, bool enabled = true) {
    std::FILE* file = enabled ? std::fopen(path.c_str(), "w") : nullptr;
    if (enabled && !file) {
        std::perror(path.c_str());
        std::exit(EXIT_FAILURE);
    }
    return OutputFile(file, CloseOutputFile{path});
}
