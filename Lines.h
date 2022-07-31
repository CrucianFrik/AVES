#ifndef AVES_LINES_H
#define AVES_LINES_H

#include "libs.h"
#include "structs.h"
#include "funcs.h"
extern Log logs;

//base calass for trajectory construction classes
class Line {
protected:
    double path_step; // длина шага цепи в метрах
    std::string files_addres;
    std::vector<Vec3D> cruve;
    std::vector<Vec3D> points;

public:
    std::vector<Vec3D> get_path() {
        return points;
    }

    virtual void build() {}
};

//calass for polyline-trajectory building
class Polyline : public Line {
public:
    Polyline(const std::vector<Vec3D> &points_, std::string f_add, double ps = 10) {
        points = points_;
        files_addres = f_add;
        path_step = ps;

        if (points_.size() < 2) {
            throw "InitError: at least 2 points are required to construct the curve";
        }
        logs.start_process("POPYLINE BUILDING");
        build();
        logs.finish_process();
    }

    Polyline() = default;

private:
    void build() override {
        for (int j = 0; j < points.size() - 1; j++) {
            int sampling = points[j].lengh(points[j + 1]) / path_step;
            if (!sampling)
                sampling = 1;

            double dx = (points[j + 1].x - points[j].x) / sampling;
            double dy = (points[j + 1].y - points[j].y) / sampling;
            double dh = (points[j + 1].h - points[j].h) / sampling;

            for (int i = 0; i <= sampling; ++i) {
                cruve.emplace_back(points[j].x + i * dx, points[j].y + i * dy, points[j].h + i * dh);
            }
        }

        std::ofstream file;
        file.open(files_addres);
        for (auto &i: cruve)
            file << std::to_string(i.x) << " " << std::to_string(i.y) << " " << std::to_string(i.h) << std::endl;
        file.close();
    }
};

//calass for CatmullROM-curve-trajectory building
class CatmullROM : public Line {
public:
    CatmullROM(const std::vector<Vec3D> &points_, std::string f_add, double min_R_, double ps = 10)
    : start_p(points_[0]), end_p(points_[points_.size() -1]), min_R(min_R_) {
        points = points_;
        files_addres = std::move(f_add);
        path_step = ps;

        if (points_.size() < 2) {
            throw "InitError: at least 2 points are required to construct the curve";
        }
        logs.start_process("CATMULLROM BUILDING");

        if (min_R)
            control_R();
        build();
        for (auto & point : points) {
            logs.add_point(point.x, point.y);
        }
        logs.finish_process();
    }

    CatmullROM() = default;

private:
    Vec3D start_p, end_p;
    double min_R;

    void build() override;

    void control_R();

    static void init_cubic_poly(double x0, double x1, double t0, double t1, CubicPoly &p) {
        p.c0 = x0;
        p.c1 = t0;
        p.c2 = -3 * x0 + 3 * x1 - 2 * t0 - t1;
        p.c3 = 2 * x0 - 2 * x1 + t0 + t1;
    }

    static void init_nonuniform_CatmullRom(double x0, double x1, double x2, double x3, double dt0, double dt1, double dt2,
                                         CubicPoly &p);

    static double vec_dist_squared(const Vec3D &p, const Vec3D &q) {
        double dx = q.x - p.x;
        double dy = q.y - p.y;
        return dx * dx + dy * dy;
    }

    static void
    init_centripetal_CR(const Vec3D &p0, const Vec3D &p1, const Vec3D &p2, const Vec3D &p3, CubicPoly &px, CubicPoly &py);

    static void get_bis_coef_and_circ_center(std::vector<Vec3D>& points, int i, double min_R_in_d_lat, double& bis_k, double& bis_b, double& x0, double& y0);

    static void get_touch_points(std::vector<Vec3D>& points, int j, double x0, double y0, double min_R_in_d_lat, Vec2D& touch_point1, Vec2D& touch_point2);

    static void check_order(std::vector<Vec3D>& points_full);

    void add_first_touch_point(std::vector<Vec3D>& points_full, Vec2D touch_point1, Vec2D touch_point2, int i);

    void add_second_touch_point(std::vector<Vec3D>& points_full, Vec2D touch_point1, Vec2D touch_point2, int i, double bis_k, double bis_b);
};

#endif //AVES_LINES_H
