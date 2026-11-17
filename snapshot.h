#pragma once

// Include after the selected variant's layer.h.
#include "run_output.h"
#ifndef MACGREGOR_REFERENCE
#include "membrane_noise.h"
#endif
#include <array>
#include <fstream>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>

class NetworkSnapshot {
    std::fstream file;
    bool reading;
    std::string path;

    void require(bool condition, const char* problem) const {
        if (!condition) throw std::runtime_error(path + ": " + problem);
    }

    void bytes(char* data, std::size_t n) {
        if (reading) file.read(data, n);
        else file.write(data, n);
        require(bool(file), reading ? "incomplete snapshot" : "snapshot write failed");
    }

    template<class T>
    void value(T& x) { bytes(reinterpret_cast<char*>(&x), sizeof(x)); }

    void signature(std::uint32_t expected) {
        auto actual = expected;
        value(actual);
        require(actual == expected, "incompatible byte order, model, topology or grid size");
    }

    void neurons(int rows, int columns, Neuron** grid) {
        for (int row = 0; row < rows; ++row) {
            for (int column = 0; column < columns; ++column) {
                auto& n = grid[row][column];
                value(n.E); value(n.Th); value(n.Gk); value(n.Ge); value(n.Gi);
                value(n.glutamate); value(n.gaba); value(n.S); value(n.P);
                require(n.S == 0 || n.S == 1, "invalid firing state");
            }
        }
    }

    // Each PSP history is mostly zero. Keep its exact length, slot positions and double values.
    void history(std::vector<double>& psps) {
        std::uint32_t length = psps.size(), used = 0;
        if (!reading)
            for (double psp : psps) if (psp != 0.0) ++used;
        value(length); value(used);
        require(length < 150 && used <= length, "invalid PSP history length");
        if (reading) {
            psps.assign(length, 0.0);
            std::uint32_t next = 0;
            for (std::uint32_t i = 0; i < used; ++i) {
                std::uint32_t slot = 0;
                value(slot);
                require(slot >= next && slot < length, "invalid PSP slot");
                value(psps[slot]);
                next = slot + 1;
            }
        } else {
            for (std::uint32_t slot = 0; slot < length; ++slot) {
                if (psps[slot] == 0.0) continue;
                value(slot); value(psps[slot]);
            }
        }
    }

    void synapses(int rows, int columns, std::vector<Synapse>** incoming, Neuron** grid) {
        for (int row = 0; row < rows; ++row) {
            for (int column = 0; column < columns; ++column) {
                auto& edges = incoming[row][column];
                std::uint32_t count = edges.size();
                value(count);
                require(count <= static_cast<unsigned>(rows * columns), "invalid incoming edge count");
                if (reading) edges.resize(count);
                for (auto& e : edges) {
                    value(e.row); value(e.clm);
                    require(e.row >= 0 && e.row < rows && e.clm >= 0 && e.clm < columns, "invalid source neuron");
                    if (reading) e.excitatory = grid[e.row][e.clm].excitatory;
                    value(e.weight); value(e.NRmax); value(e.NR); value(e.tau);
                    history(e.psps);
                }
            }
        }
    }

    void randomState() {
#ifdef WIN32
        throw std::runtime_error("Restart snapshots require the POSIX drand48 generator");
#else
        // seed48 exposes the previous 48-bit state. Copy it before restoring it.
        std::array<unsigned short, 3> seed{};
        if (!reading) {
            const auto old = seed48(seed.data());
            std::copy(old, old + 3, seed.begin());
            seed48(seed.data());
        }
        bytes(reinterpret_cast<char*>(seed.data()), seed.size() * sizeof(seed[0]));
        if (reading) seed48(seed.data());
        int workers = 0;
        std::vector<std::string> streams;
        if (!reading) {
#pragma omp parallel
            {
#pragma omp single
                { workers = omp_get_num_threads(); streams.resize(workers); }
#ifndef MACGREGOR_REFERENCE
                // The reference recreates its normal generator on every call.
                auto& noise = membraneNoise();
                std::ostringstream out;
                out.imbue(std::locale::classic());
                out << noise.generator << ' ' << noise.distribution;
                streams[omp_get_thread_num()] = out.str();
#endif
            }
        }
        value(workers);
        require(workers > 0 && workers <= omp_get_thread_limit() && workers <= 65536, "invalid OpenMP team size");
        if (reading) {
            omp_set_dynamic(0);
            omp_set_num_threads(workers);
        }
#ifndef MACGREGOR_REFERENCE
        if (reading) streams.resize(workers);
        for (auto& stream : streams) {
            std::uint32_t n = stream.size();
            value(n);
            require(n && n < 100000, "invalid normal generator state size");
            if (reading) stream.resize(n);
            bytes(&stream[0], n);
        }
        if (reading) {
            int valid = 1;
#pragma omp parallel reduction(&:valid)
            {
                std::istringstream in(streams[omp_get_thread_num()]);
                in.imbue(std::locale::classic());
                auto& noise = membraneNoise();
                valid &= omp_get_num_threads() == workers && bool(in >> noise.generator >> noise.distribution);
            }
            require(valid != 0, "cannot restore OpenMP normal generators");
        }
#endif
#endif
    }

public:
    NetworkSnapshot(const std::string& path, bool reading)
        : file(path, std::ios::binary | (reading ? std::ios::in : std::ios::out | std::ios::trunc)),
          reading(reading), path(path) {
        require(bool(file), "cannot open snapshot");
        static_assert(sizeof(int) == 4 && sizeof(unsigned short) == 2 && sizeof(double) == 8 &&
                      std::numeric_limits<double>::is_iec559, "unsupported snapshot scalar format");
    }

    void header(Topology topology, bool reference, int rows, int columns, int& nextStep, Param& p) {
        char magic[] = "MGSNAP01";
        bytes(magic, 8);
        require(std::memcmp(magic, "MGSNAP01", 8) == 0, "unknown snapshot format");
        signature(0x01020304); signature(static_cast<unsigned>(topology)); signature(reference);
        signature(rows); signature(columns); value(nextStep);
        require(nextStep > 0, "invalid next step");
        value(p.ExTgk); value(p.ExB); value(p.ExTth); value(p.ExTe1); value(p.ExTe2); value(p.ExTi1); value(p.ExTi2);
        value(p.InTgk); value(p.InB); value(p.InTth); value(p.InTe1); value(p.InTe2); value(p.InTi1); value(p.InTi2);
        value(p.TauRelease); value(p.TauReplenish); value(p.SpontRelease); value(p.Neighborhood); value(p.NRmax);
    }

    void state(int rows, int columns, Neuron** grid, std::vector<Synapse>** incoming,
               int nextStep, EpochState& epoch) {
        neurons(rows, columns, grid);
        synapses(rows, columns, incoming, grid);
        randomState();
        value(epoch.interval); value(epoch.firstStep);
        require(epoch.interval >= 0 && epoch.firstStep >= 0 && epoch.firstStep <= nextStep, "invalid epoch");
        if (epoch.interval) {
            require(nextStep - epoch.firstStep < epoch.interval, "invalid partial epoch length");
            value(epoch.excitatory); value(epoch.inhibitory); value(epoch.fields.legacy); value(epoch.fields.fixed);
            if (reading) { epoch.spikes.resize(rows * columns); epoch.glutamate.resize(rows * columns); }
            bytes(reinterpret_cast<char*>(epoch.spikes.data()), epoch.spikes.size() * sizeof(int));
            bytes(reinterpret_cast<char*>(epoch.glutamate.data()), epoch.glutamate.size() * sizeof(double));
        }
        if (reading) require(file.peek() == EOF && !file.bad(), "unexpected data after snapshot");
        else {
            file.close();
            require(bool(file), "could not finish snapshot");
        }
    }
};
