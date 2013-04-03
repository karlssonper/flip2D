#ifndef LOG_H_
#define LOG_H_

#include <iostream>
#include <fstream>
#include <sstream>

inline std::string NowTime();

#define LOG_OUTPUT(args)                        \
    Log::instance().stream() << "-- " <<        \
        NowTime() << ": " << args;              \
    Log::instance().output();                   \

#define LOG_OUTPUT_WITHOUT_TIMESTAMPS(args)     \
    Log::instance().stream() << args;           \
    Log::instance().output();                   \

#define LOG_ERROR(args)                         \
    Log::instance().stream() << "-- " <<        \
        NowTime() << ": ERROR-- " <<args;       \
    Log::instance().error();                    \
    
#define LOG_DEBUG(args)                         \
    Log::instance().stream() << args;           \
    Log::instance().output();                   \

class Log
{
  public:
    void error()
    {
        std::cout << _ss.str() << std::endl;
        _logfile << _ss.str() << std::endl;
        _ss.str("");
        _ss.clear();
    }

    void output()
    {
        std::cerr << _ss.str() << std::endl;
        _logfile << _ss.str() << std::endl;
        _ss.str("");
        _ss.clear();
    }

    static Log & instance()
    {
        static Log instance;
        return instance;
    }

    std::stringstream & stream()
    {
        return _ss;
    }
    
  protected:
    std::stringstream _ss;
    std::ofstream _logfile;

    Log() : _logfile("flip2D_sim.log")
    {
        if (!_logfile.is_open()) {
            std::cerr << "Could not create log file" << std::endl;
        }
    }
    
    Log(const Log &);
    void operator=(const Log &);
};

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__)

#include <windows.h>

inline std::string NowTime()
{
    const int MAX_LEN = 200;
    char buffer[MAX_LEN];
    if (GetTimeFormatA(LOCALE_USER_DEFAULT, 0, 0, 
            "HH':'mm':'ss", buffer, MAX_LEN) == 0)
        return "Error in NowTime()";

    char result[100] = {0};
    static DWORD first = GetTickCount();
    std::sprintf(result, "%s.%03ld", buffer, (long)(GetTickCount() - first) % 1000); 
    return result;
}

#else

#include <sys/time.h>

inline std::string NowTime()
{
    char buffer[11];
    time_t t;
    time(&t);
    tm r = {0};
    strftime(buffer, sizeof(buffer), "%X", localtime_r(&t, &r));
    struct timeval tv;
    gettimeofday(&tv, 0);
    char result[100] = {0};
    std::sprintf(result, "%s.%03ld", buffer, (long)tv.tv_usec / 1000); 
    return result;
}

#endif //WIN32

#endif
