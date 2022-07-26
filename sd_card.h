#ifndef SD_H
#define SD_H

#include <EEPROM.h>
#define EEPROM_SIZE 1

#include "FS.h"
#include "SD.h"
#include "SPI.h"

#define MISO  14
#define SCK  27
#define MOSI  26
#define CS  25

int filenumber;
String filestr;
const char* filename;

SPIClass spi = SPIClass(VSPI);

const char HEADER[] = {
    "ch1,"
    "ch2,"
    "ch3,"
    "ch4,"
    "ch5,"
    "ch6,"
    "ch7,"
    "ch8,"
    "ch9,"
    "ch10,"
    "roll,"
    "pitch,"
    "yaw,"
    "Altitude,"
    "Latitude,Longitude,"
    "time"
};


bool writeFile(fs::FS &fs, const char * path, const char * message){
  File file = fs.open(path, FILE_WRITE);
  if(!file){
    return false;
  }
  file.print(message);
  file.close();
  return true;
}

void appendFile(fs::FS &fs, const char * path, const char * message){
  File file = fs.open(path, FILE_APPEND);
  if(!file){
    return;
  }
  file.print(message);
  file.close();
}

bool init_sd(){
  EEPROM.begin(EEPROM_SIZE);
  filenumber = EEPROM.read(0);
  if ((filenumber<0) or (filenumber>256)){
    filenumber=0;
  }
  filestr="/data_"+String(filenumber)+".txt";
  filename=filestr.c_str();
  spi.begin(SCK, MISO, MOSI, CS);
  SD.begin(CS,spi,80000000);
  EEPROM.write(0, filenumber+1);
  if(!writeFile(SD, filename, HEADER)){
    return false;
  }
  EEPROM.commit();
  return true;
}

void to_file(String message){
  appendFile(SD, filename, message.c_str());
  }
#endif
