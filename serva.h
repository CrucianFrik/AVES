#ifndef SERVA_H
#define SERVA_H

#include <ESP32Servo.h>
#include "FlySkyIBus.h"

const int aileron_l_Pin = 18;
const int aileron_r_Pin = 23;
const int pgo_l_pin = 5;
const int pgo_r_pin = 4;
const int antenaPin  = 35;
const int engine_l_Pin  = 13;
const int engine_r_Pin  = 19;



HardwareSerial FSkySerial(1);

//Petr
double roll_ctrl_effect, pitch_ctrl_effect;
bool hand_control_flag = true;
/* servo_1 - pgo left
   servo_2 - pgo right
   servo_3, servo_4 - engines 
   servo_5 - aileron_l
   servo_6 - aileron_r*/
/// end Petr///

Servo servo_1;  
Servo servo_2;  
Servo servo_3;
Servo servo_4;
Servo servo_5;
Servo servo_6;

int servo_control[10] = {1500,1500,1000,1000,1500,1500,1500,1500,1500,1500};
int servo_control_buf;

void init_control() {
  IBus.begin(FSkySerial);
  servo_1.attach(pgo_l_pin); 
  servo_2.attach(pgo_r_pin); 
  servo_3.attach(engine_l_Pin); 
  servo_4.attach(engine_r_Pin); 
  servo_5.attach(aileron_l_Pin); 
  servo_6.attach(aileron_r_Pin);
}

void read_chanel(int channel){
  int servo_control_buf = IBus.readChannel(channel);
//  if (servo_control_buf>=1000 && servo_control_buf<=2000)
  {
    servo_control[channel] = servo_control_buf;
  }
}

void read_control() {
  IBus.loop();
  read_chanel(0);
  read_chanel(1);
  read_chanel(2);
  read_chanel(3);
  read_chanel(4);
  read_chanel(5);
  read_chanel(6);
  read_chanel(7);
  read_chanel(8);
  read_chanel(9);
  hand_control_flag = (servo_control[6] > 1500)? true : false;
}

void write_control_hand_petr(){
   servo_1.write(servo_control[1]);
   servo_2.write(3000-servo_control[1]);
   servo_3.write(servo_control[2]);
   servo_4.write(servo_control[3]);
   servo_5.write(servo_control[0]);
   servo_6.write(servo_control[0]);
}

//void write_control_pid_petr(){
//  servo_1.write(1000);
//}

#endif
