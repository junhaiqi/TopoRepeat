#ifndef LOG_H
#define LOG_H

#include <chrono>
#include <iomanip>
#include <iostream>
#include <sys/resource.h>
#include <sys/time.h>

// Global timer start point
static std::chrono::steady_clock::time_point _global_start_time;

// Call this function at the very beginning of main()
inline void start_main_timer()
{
    _global_start_time = std::chrono::steady_clock::now();
}

// Helper function: Get current time and print formatted resource usage log
// level: Log level tag (e.g., "INFO", "main", "DONE")
inline void print_resource_usage(const std::string &level = "INFO")
{
    // --- 1. Calculate Resource Usage ---

    // Calculate Real Time (Wall-clock time)
    auto end_time = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff = end_time - _global_start_time;
    double real_time = diff.count();

    // Get system resource usage (CPU & Memory)
    struct rusage r;
    getrusage(RUSAGE_SELF, &r);

    // Calculate CPU time (User time + System time)
    double cpu_time = (r.ru_utime.tv_sec + r.ru_utime.tv_usec / 1.0e6) +
                      (r.ru_stime.tv_sec + r.ru_stime.tv_usec / 1.0e6);

    // Calculate Peak RSS (Resident Set Size)
    // Note: r.ru_maxrss is in Kilobytes on Linux, but Bytes on macOS.
    double peak_rss_gb = 0.0;
#ifdef __APPLE__
    peak_rss_gb = (double)r.ru_maxrss / 1024.0 / 1024.0 / 1024.0;
#else
    // Linux default (KB to GB)
    peak_rss_gb = (double)r.ru_maxrss / 1024.0 / 1024.0;
#endif

    // --- 2. Build the Message String ---
    std::stringstream ss;
    ss << std::fixed << std::setprecision(3);
    ss << "Real time: " << real_time << " sec; "
       << "CPU: " << cpu_time << " sec; "
       << "Peak RSS: " << peak_rss_gb << " GB";

    std::string message = ss.str();

    // --- 3. Get Current Calendar Time (for Log Timestamp) ---
    auto now_log = std::chrono::system_clock::now();
    std::time_t log_time_t = std::chrono::system_clock::to_time_t(now_log);

    // --- 4. Print in the Specified Format ---
    // Format: [YYYY-MM-DD HH:MM:SS] [LEVEL] Message
    // std::cerr << "[" << std::put_time(std::localtime(&log_time_t), "%Y-%m-%d %H:%M:%S") << "] "
    //           << "[" << level << "] "
    //           << message << std::endl;
    std::cerr << "[" << level << "] "
              << message << std::endl;
}

inline void log(const std::string &message, const std::string &level = "INFO")
{
    auto now = std::chrono::system_clock::now();
    auto time_t = std::chrono::system_clock::to_time_t(now);

    std::cerr << "[" << std::put_time(std::localtime(&time_t), "%Y-%m-%d %H:%M:%S") << "] "
              // << "[" << std::setw(5) << level << "] "
              << "[" << level << "] "
              << message << std::endl;
}

#endif