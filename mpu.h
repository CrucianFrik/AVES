#ifndef MPU_H
#define MPU_H


#include "I2Cdev.h"
#include "MPU9250_9Axis_MotionApps41.h"
#include "Wire.h"

MPU9250 mpu;

uint8_t mpuIntStatus,devStatus,fifoBuffer[64];
uint16_t packetSize,fifoCount,mag[3];
float f_mag[3],d_mag[4];
Quaternion q,q_mag;
VectorInt16 aa, gg;

double roll,pitch,yaw;

bool init_mpu() {
//    Wire.begin();
//    Wire.setClock(400000);
    mpu.initialize();
    devStatus = mpu.dmpInitialize();
    mpu.setDMPEnabled(true);
    mpuIntStatus = mpu.getIntStatus();
    packetSize = mpu.dmpGetFIFOPacketSize();
    return true;
}

void update_mpu() {
    mpuIntStatus = mpu.getIntStatus();
    fifoCount = mpu.getFIFOCount();
    if ((mpuIntStatus & 0x10) || fifoCount == 1024) {
        mpu.resetFIFO();
    } else if (mpuIntStatus & 0x02) {
        while (fifoCount < packetSize) fifoCount = mpu.getFIFOCount();
        mpu.getFIFOBytes(fifoBuffer, packetSize);
        fifoCount -= packetSize;
        mpu.dmpGetQuaternion(&q, fifoBuffer);
        mpu.dmpGetAccel(&aa, fifoBuffer);
        mpu.dmpGetGyro(&gg, fifoBuffer);
        roll=atan2(2*(q.w*q.x+q.y*q.z),1-2*(q.x*q.x+q.y*q.y))* 180/M_PI;
        pitch=-asin(2*(q.w*q.y-q.x*q.z))* 180/M_PI;
        yaw=atan2(2*(q.w*q.z+q.x*q.y),1-2*(q.y*q.y+q.z*q.z))* 180/M_PI;
    }
}

#endif
