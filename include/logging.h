#pragma once

#include <mpi.h>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_sinks.h>

#include <stdexcept>
#include <string>

namespace Log {

struct Config {
    spdlog::level::level_enum level = spdlog::level::info;
};

inline Config parse_options(int argc, char *argv[])
{
    Config config;
    for (int i = 3; i < argc; ++i) {
        const std::string option = argv[i];
        if (option != "--log-level")
            throw std::invalid_argument("Unknown option '" + option + "'");
        if (++i == argc)
            throw std::invalid_argument(option + " requires a value");
        config.level = spdlog::level::from_str(argv[i]);
    }
    return config;
}

inline void init(const Config &config)
{
    if (spdlog::get("FANS"))
        return;

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    auto logger = spdlog::stdout_logger_mt("FANS");
    logger->set_pattern(config.level <= spdlog::level::debug
                            ? "[%Y-%m-%d %H:%M:%S.%e] [%n] [%l] [pid %P] %v"
                            : "%v");
    logger->set_level(rank == 0 ? config.level : spdlog::level::off);
}

inline spdlog::logger &logger()
{
    return *spdlog::get("FANS");
}

} // namespace Log
