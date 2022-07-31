#include "GPS_Path.h"
#include <iostream>

Log logs{};

int main() {
    //logs.init("" + add_time("LOG") + ".txt");
    GPS_Path GPS{"", "GPS_text.txt", PathType::curve, 5, 200, 5, 20};
}
