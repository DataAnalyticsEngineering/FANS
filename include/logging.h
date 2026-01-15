#ifndef LOGGING_H
#define LOGGING_H

#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <memory>
#include "mpi.h"

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

class Logger;
class MPI_TraceSync {
  public:
    explicit MPI_TraceSync(Logger &log, bool append);
    ~MPI_TraceSync();

    template <class T>
    MPI_TraceSync &operator<<(const T &v);
    MPI_TraceSync &operator<<(std::ostream &(*m)(std::ostream &) );

  private:
    Logger            &_log;
    bool               _append;
    std::ostringstream _buffer;
};

template <class T>
MPI_TraceSync &MPI_TraceSync::operator<<(const T &v)
{
    _buffer << v;
    return *this;
}

class Logger {
  public:
    friend class MPI_TraceSync;
    explicit Logger(std::string prefix, int comm_rank, int comm_size, const MPI_Comm &comm);

    /// error > info > warn > debug > trace
    std::ostream &error(bool append = false);
    /// error > info > warn > debug > trace
    std::ostream &info(bool append = false);
    /// error > info > warn > debug > trace
    std::ostream &warn(bool append = false);
    /// error > info > warn > debug > trace
    std::ostream &debug(bool append = false);
    /// error > info > warn > debug > trace
    MPI_TraceSync trace(bool append = false);
    /// progress bar
    void progress(const std::string &prefix, int step, int max);

  private:
    std::ostream &trace_impl(bool append = false);
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
    /// communicator
    MPI_Comm _comm;
};

/**
 * Creates all loggers and sets the level
 * */
void init(int comm_rank, int comm_size, const MPI_Comm &comm);

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
