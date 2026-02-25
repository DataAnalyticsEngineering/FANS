#include "logging.h"

#include <utility>

#include "mpi.h"

int active_rank = 0;

Log::MPI_TraceSync::MPI_TraceSync(Logger &log, bool append)
    : _log(log), _append(append) {}
Log::MPI_TraceSync::~MPI_TraceSync()
{
    if (active_level >= All) {
        MPI_Barrier(_log._comm);
        const int size = _log._comm_size;
        for (int i = 0; i < size; i++) {
            Log::setActiveRank(i);
            MPI_Barrier(_log._comm);
            _log.trace_impl(_append) << _buffer.str() << std::flush;
        }
        Log::setActiveRank(0);
        MPI_Barrier(_log._comm);
    }
}

Log::MPI_TraceSync &Log::MPI_TraceSync::operator<<(std::ostream &(*m)(std::ostream &) )
{
    _buffer << m;
    return *this;
}

Log::Logger::Logger(std::string prefix, int comm_rank, int comm_size, const MPI_Comm &comm)
    : _prefix(std::move(prefix)), _nullstr(), _comm_rank(comm_rank), _comm_size(comm_size), _comm(comm)
{
    _start_time = std::chrono::steady_clock::now();
    _nullstr.setstate(std::ios_base::badbit);
}

std::ostream &Log::Logger::error(bool append)
{
    if (active_level >= Error && _comm_rank == active_rank) {
        if (append)
            return std::cout;
        std::cout << _prefix;
        auto now         = std::chrono::steady_clock::now();
        auto elapse_time = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(now - _start_time).count()) / 1000.0;
        std::cout << elapse_time << ": ";
        return std::cout;
    } else
        return _nullstr;
}

std::ostream &Log::Logger::info(bool append)
{
    if (active_level >= Info && _comm_rank == active_rank) {
        if (append)
            return std::cout;
        std::cout << _prefix;
        auto now         = std::chrono::steady_clock::now();
        auto elapse_time = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(now - _start_time).count()) / 1000.0;
        std::cout << elapse_time << ": ";
        return std::cout;
    } else
        return _nullstr;
}

std::ostream &Log::Logger::warn(bool append)
{
    if (active_level >= Warn && _comm_rank == active_rank) {
        if (append)
            return std::cout;
        std::cout << _prefix;
        auto now         = std::chrono::steady_clock::now();
        auto elapse_time = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(now - _start_time).count()) / 1000.0;
        std::cout << elapse_time << ": ";
        return std::cout;
    } else
        return _nullstr;
}

std::ostream &Log::Logger::debug(bool append)
{
    if (active_level >= Debug && _comm_rank == active_rank) {
        if (append)
            return std::cout;
        std::cout << _prefix;
        auto now         = std::chrono::steady_clock::now();
        auto elapse_time = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(now - _start_time).count()) / 1000.0;
        std::cout << elapse_time << ": ";
        return std::cout;
    } else
        return _nullstr;
}

Log::MPI_TraceSync Log::Logger::trace(bool append)
{
    return Log::MPI_TraceSync(*this, append);
}

std::ostream &Log::Logger::trace_impl(bool append)
{
    if (active_level >= All && _comm_rank == active_rank) {
        if (append)
            return std::cout;
        std::cout << _prefix << "[" << active_rank << "] ";
        auto now         = std::chrono::steady_clock::now();
        auto elapse_time = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(now - _start_time).count()) / 1000.0;
        std::cout << elapse_time << ": ";
        return std::cout;
    } else
        return _nullstr;
}

void Log::Logger::progress(const std::string &prefix, int step, int max)
{
    if (_comm_rank != active_rank)
        return;

    if (step == 0)
        std::cout << _prefix << prefix;

    int digits = 1;
    int div    = 10;
    while ((max / div) > 0) {
        digits++;
        div *= 10;
    }
    if (step > 0)
        for (int i = 0; i < digits; i++)
            std::cout << '\b';

    int step_digits = 1;
    div             = 10;
    while ((step / div) > 0) {
        step_digits++;
        div *= 10;
    }

    for (int i = 0; i < digits - step_digits; i++)
        std::cout << ' ';
    std::cout << step;

    if (step >= max - 1) {
        for (int i = 0; i < digits; i++)
            std::cout << '\b';
        for (int i = 0; i < prefix.size(); i++)
            std::cout << '\b';
        for (int i = 0; i < _prefix.size(); i++)
            std::cout << '\b';
    }
    // if (step >= max) std::cout << '\n';
    std::cout.flush();
}

void Log::init(int comm_rank, int comm_size, const MPI_Comm &comm)
{
    Level level = Log::Error;
    if constexpr (VERBOSITY <= 0)
        level = Error;
    if constexpr (VERBOSITY == 1)
        level = Info;
    if constexpr (VERBOSITY == 2)
        level = Warn;
    if constexpr (VERBOSITY == 3)
        level = Debug;
    if constexpr (VERBOSITY >= 4)
        level = All;

    active_level = level;
    general      = std::make_unique<Logger>("[FANS-GENERAL] ", comm_rank, comm_size, comm);
    solver       = std::make_unique<Logger>("[FANS-SOLVER] ", comm_rank, comm_size, comm);
    io           = std::make_unique<Logger>("[FANS-IO] ", comm_rank, comm_size, comm);
}

void Log::finalize()
{
    general = nullptr;
    solver  = nullptr;
    io      = nullptr;
}

Log::Level                   Log::active_level = Info;
std::unique_ptr<Log::Logger> Log::general      = nullptr;
std::unique_ptr<Log::Logger> Log::solver       = nullptr;
std::unique_ptr<Log::Logger> Log::io           = nullptr;

const std::string Log::level_to_string(Log::Level level)
{
    switch (level) {
    case Error:
        return "Error";
    case Info:
        return "Info";
    case Warn:
        return "Warn";
    case Debug:
        return "Debug";
    case All:
        return "All";
    }
    return "";
}

void Log::setActiveRank(int rank)
{
    active_rank = rank;
}
