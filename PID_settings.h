#ifndef PID_SETTINGS_H
#define PID_SETTINGS_H

#include "PID.h"

#define ROLL_P 5.0
#define ROLL_I 0.0
#define ROLL_D 0.0
#define ROLL_MIN 1000.0
#define ROLL_MAX 2000.0
#define PITCH_P 5.0
#define PITCH_I 0.0
#define PITCH_D 0.0
#define PITCH_MIN 1000.0
#define PITCH_MAX 2000.0

// 2 ch 1500 1146 1854 (center, max, min)
// 5 ch 1349 989 1704
PID roll_pid(ROLL_P, ROLL_I, ROLL_D);
PID pitch_pid(PITCH_P, PITCH_I, PITCH_D);

#endif
