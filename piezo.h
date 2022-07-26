#ifndef PIEZO_H
#define PIEZO_H

int port = 12;

void init_piezo(){
  pinMode(port,OUTPUT);
}

void delay_piezo(int time_delay){
  digitalWrite(port,HIGH);
  delay(time_delay);
  digitalWrite(port,LOW);
}

#endif
