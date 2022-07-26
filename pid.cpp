#include "pid.h"
// PID::PID(double p, double i, double d, double maxOut, double minOut){
PID::PID(double p, double i, double d){

//    t = get_time();
    set_pid(p,i,d);
}

void PID::set_pid(double pp, double ii, double dd) {
    p = pp;
    i = ii;
    d = dd;
}

double PID::ctrl(double req, double cur) {
    double res;
    e = req - cur;
//    sum += (get_time()() - t)*e;
    res = p*e; // + d*(e - prev)/t + sum*i;
    prev  = e;
    // return ((res > maxOut)? maxOut : (res < minOut)? minOut : res);
    return res;
}
void PID::reset(){
    e = 0;
    prev = 0;
    sum = 0;
}
