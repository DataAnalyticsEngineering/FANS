#pragma once

#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <memory>
#include <vector>
#include "mpi.h"

#include <sys/stat.h>

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

    /// error > warn > debug > info > trace
    std::ostream &error(bool append = false);
    /// error > warn > debug > info > trace
    std::ostream &info(bool append = false);
    /// error > warn > debug > info > trace
    std::ostream &warn(bool append = false);
    /// error > warn > debug > info > trace
    std::ostream &debug(bool append = false);
    /// error > warn > debug > info > trace
    MPI_TraceSync trace(bool append = false);
    /// progress bar
    void progress(const std::string &prefix, int step, int max) const;

  private:
    std::ostream &trace_impl(bool append = false);
    /// starting time
    std::chrono::steady_clock::time_point _start_time;
    /// computes the elapsed time in (hours, minutes, seconds)
    std::tuple<int, int, int> get_elapsed_time() const;
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
 * Formats a string according to fmt with the provided args
 * @param fmt format string
 * @param args format string arguments
 * @return formatted string
 */
template <typename... Args>
std::string format(const std::string &fmt, Args &&...args)
{
    const int         buf_size = std::snprintf(nullptr, 0, fmt.c_str(), args...);
    std::vector<char> buf(buf_size + 1);
    std::snprintf(buf.data(), buf.size(), fmt.c_str(), args...);
    return {buf.data(), static_cast<std::string::size_type>(buf_size)};
}

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
