#ifndef LOGGING_H
#define LOGGING_H

#include <chrono>
#include <fstream>
#include <iostream>
#include <memory>

namespace Log {

enum Level {
    Error,
    Info,
    Warn,
    Debug,
    All
};

[[maybe_unused]] const std::string level_to_string(Level level);

extern Level active_level;

class Logger {
  public:
    explicit Logger(std::string prefix, int comm_rank, int comm_size);

    /// error > info > warn > debug > trace
    std::ostream& error(bool append=false);
    /// error > info > warn > debug > trace
    std::ostream& info(bool append=false);
    /// error > info > warn > debug > trace
    std::ostream& warn(bool append=false);
    /// error > info > warn > debug > trace
    std::ostream& debug(bool append=false);
    /// error > info > warn > debug > trace
    std::ostream& trace(bool append=false);
    /// progress bar
    void progress(const std::string& prefix, int step, int max);

  private:
    /// starting time
    std::chrono::steady_clock::time_point _start_time;
    /// what the logger should always print first
    const std::string _prefix;
    /// empty stream to write nothing
    std::ofstream _nullstr;
    /// MPI comm rank
    int _comm_rank;
    /// MPI comm size
    int _comm_size;

};

/**
 * Creates all loggers and sets the level
 * */
void init(int comm_rank, int comm_size);

/**
 * Frees all memory
 * */
void finalize();

/**
 * Set activate rank
 */
[[maybe_unused]] void setActiveRank(int rank);

/// logger with prefix [GENERAL]
extern std::unique_ptr<Logger> general;
/// logger with prefix [SOLVER]
extern std::unique_ptr<Logger> solver;
/// logger with prefix [IO]
extern std::unique_ptr<Logger> io;

}; // namespace Log

#endif