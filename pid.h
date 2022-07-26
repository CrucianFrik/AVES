#ifndef PID_H
#define PID_H

class PID {
private:
    double prev = 0; //previous error value
    double p = 0;    //p, i, d - coefficients
    double i = 0;
    double d = 0;
    double sum = 0;    //integral sum
    double e = 0;    //error (difference of current and required values)
    double t;       //time step
    double maxOut;
    double minOut;

public:
    // PID(double p, double i, double d, double maxOut, double minOut);     //it is necessary to define time step
    PID(double p, double i, double d); 
    double ctrl(double req, double cur); // req - required; cur - current;
    void set_pid(double pp, double ii, double dd);
    void reset();
    void set_borders();
};

#endif