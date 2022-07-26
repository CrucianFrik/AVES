#ifndef BMP180_H
#define BMP180_H
#include <Wire.h>

#include <Adafruit_BMP085.h>
Adafruit_BMP085 bmp;

bool init_bmp(){
  #ifdef I2C
  #define I2C
  Wire.begin();
  #endif
  return bmp.begin(0);
}

float alt()
{
  return bmp.readAltitude();
}

#endif
