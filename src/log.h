#ifndef LOG_H_
#define LOG_H_

#include <iostream>
#inclide <fstream>

class Log
{
    std::ostringstream error()
    {
        os << "- " << NowTime();
        os << " " << "FLIP2D ERROR: ";
        return os;
    }
    
  protected:
    static Log log;

    std::ostringstream os;
    
    Log();
    Log(const Log &);
    void operator=(const Log &);
};

#endif
