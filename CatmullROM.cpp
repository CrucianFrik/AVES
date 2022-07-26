#include "Lines.h"

void CatmullROM::build() {
    for (int j = 0; j < points.size() - 1; j++) {
        int sampling = points[j].lengh(points[j + 1]) / path_step;
        if (!sampling)
            sampling = 1;
        auto step = double(1. / sampling);

        CubicPoly px{}, py{};
        CubicPoly pi{}, ph{};

        Vec3D p0, p3;
        if (!j)
            p0 = start_p;
        else
            p0 = points[j - 1];

        if (j == points.size() - 2)
            p3 = end_p;
        else
            p3 = points[j + 2];

        auto jj = static_cast<double>(j);
        init_centripetal_CR(Vec3D{jj, p0.h, 1}, Vec3D{jj + 0.1, points[j].h, 0}, Vec3D{jj + 0.2, points[j + 1].h, 0},
                          Vec3D{jj + 0.3, p3.h, 0}, pi, ph);
        init_centripetal_CR(p0, points[j], points[j + 1], p3, px, py);
        for (int i = 0; i <= sampling; ++i)
            cruve.emplace_back(px.eval(step * i), py.eval(step * i), ph.eval(step * i));
    }

    std::ofstream file;
    file.open(files_addres);
    for (auto &i: cruve)
        file << std::to_string(i.x) << " " << std::to_string(i.y) << " " << std::to_string(i.h) << std::endl;
    file.close();
}

void CatmullROM::init_centripetal_CR(const Vec3D &p0, const Vec3D &p1, const Vec3D &p2, const Vec3D &p3,
                                   CubicPoly &px, CubicPoly &py) {
    double dt0 = powf(vec_dist_squared(p0, p1), 0.25f);
    double dt1 = powf(vec_dist_squared(p1, p2), 0.25f);
    double dt2 = powf(vec_dist_squared(p2, p3), 0.25f);

    // safety check for repeated points
    if (dt1 < 1e-4f) dt1 = 1.0f;
    if (dt0 < 1e-4f) dt0 = dt1;
    if (dt2 < 1e-4f) dt2 = dt1;

    init_nonuniform_CatmullRom(p0.x, p1.x, p2.x, p3.x, dt0, dt1, dt2, px);
    init_nonuniform_CatmullRom(p0.y, p1.y, p2.y, p3.y, dt0, dt1, dt2, py);
}

void
CatmullROM::init_nonuniform_CatmullRom(double x0, double x1, double x2, double x3, double dt0, double dt1, double dt2,
                                     CubicPoly &p) {
// compute coefficients for a nonuniform Catmull-Rom spline
// compute tangents when parameterized in [t1,t2]
    double t1 = (x1 - x0) / dt0 - (x2 - x0) / (dt0 + dt1) + (x2 - x1) / dt1;
    double t2 = (x2 - x1) / dt1 - (x3 - x1) / (dt1 + dt2) + (x3 - x2) / dt2;

// rescale tangents for parametrization in [0,1]
    t1 *= dt1;
    t2 *= dt1;

    init_cubic_poly(x1, x2, t1, t2, p);
}
