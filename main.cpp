#include "path_creation/GPS_Path.h"
#include <iostream>

Log logs{};

int main() {
    logs.init("/logs/" + add_time("LOG") + ".txt");
    GPS_Path GPS{"/", "GPS_text.txt", PathType::curve, 5, 5, 10};
}
