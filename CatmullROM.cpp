#include "Lines.h"


void CatmullROM::control_R(){
    std::ofstream center_file;
    center_file.open("centers.txt");

    std::ofstream touch_point_file;
    touch_point_file.open("touch_points.txt");

    double deg_lon = 111003;
    double deg_lat = cos(points[0].x)*6360*2*M_PI*1000/360; // 1000 -- km to m
    double min_R_in_d_lat = ((min_R / deg_lon) + (min_R / deg_lat))/2;


    std::cout << deg_lat << "\n";

    if (points.size() <= 2)
        return;

    std::vector<Vec3D> points_full;
    points_full.push_back(points[0]);

    for (int i = 1; i < points.size() - 1 ; i++) {
        double rho01 = points[i].lengh(points[i - 1]);
        double rho21 = points[i].lengh(points[i + 1]);

        Vec2D V10 = Vec2D{(points[i - 1].x - points[i].x), (points[i - 1].y - points[i].y)} / rho01;
        Vec2D V12 = Vec2D{(points[i + 1].x - points[i].x), (points[i + 1].y - points[i].y)} / rho21;

        Vec2D point12_ = V10*20 + points[i].get_2d();
        Vec2D point12 = V12*20 + points[i].get_2d();
        Vec2D mid_point = point12 + (point12_ - point12) / 2;

        //уравнение биссектрисы:
        double mid_k = (points[i].y - mid_point.y) / (points[i].x - mid_point.x);
        double mid_b = points[i].y - mid_k*points[i].x;

        double rho1m = mid_point.lengh(points[i].get_2d());
        double len = sqrt(pow(mid_point.x - points[i].get_2d().x, 2) + pow(mid_point.y - points[i].get_2d().y, 2));
        Vec2D V1m = ((mid_point - points[i].get_2d()) / len) * (min_R_in_d_lat);

        Vec2D center = V1m + points[i].get_2d();

        //////////////////////////////////////
        const auto digits = std::numeric_limits<double>::digits10;
        center_file << std::setfill(' ') << std::setw(digits + 4);
        center_file << std::fixed << std::setprecision(digits) << center.x << " ";
        center_file << std::setfill(' ') << std::setw(digits + 4);
        center_file << std::fixed << std::setprecision(digits) << center.y << " ";
        center_file << "\n";
        //////////////////////////////////////

        double x0 = center.x, y0 = center.y;

        for (int j = i - 1; j <= i + 1; j+=2) {
            double x1 = points[j].x, y1 = points[j].y;
            double k_1, k_2, b_1, b_2;
            Vec2D touch_point1, touch_point2;

            double dx = x1-x0;
            double dy = y1-y0;

            if (x0 - x1 == min_R_in_d_lat) {
                touch_point1 = Vec2D{x0 - min_R_in_d_lat, y0};
                k_2 = -(pow(dy, 2) - pow(min_R_in_d_lat, 2)) / (2 * min_R_in_d_lat * dy);
                b_2 = y1 - k_2 * x1;
                double x = (2 * x0 - 2 * k_2 * b_2 + 2 * k_2 * y0) / (2 * (1 + pow(k_2, 2)));
                touch_point2 = Vec2D{x, k_2 * x + b_2};
            } else if (x0 - x1 == -min_R_in_d_lat) {
                touch_point1 = Vec2D{x0 + min_R_in_d_lat, y0};
                k_2 = (pow(dy, 2) - pow(min_R_in_d_lat, 2)) / (2 * min_R_in_d_lat * dy);
                b_2 = y1 - k_2 * x1;
                double x = (2 * x0 - 2 * k_2 * b_2 + 2 * k_2 * y0) / (2 * (1 + pow(k_2, 2)));
                touch_point2 = Vec2D{x, k_2 * x + b_2};
            } else {
                k_1 = (dx * dy + min_R_in_d_lat * sqrt(pow(dx, 2) + pow(dy, 2) - pow(min_R_in_d_lat, 2))) / (pow(dx, 2) - pow(min_R_in_d_lat, 2));
                k_2 = (dx * dy - min_R_in_d_lat * sqrt(pow(dx, 2) + pow(dy, 2) - pow(min_R_in_d_lat, 2))) / (pow(dx, 2) - pow(min_R_in_d_lat, 2));
                b_1 = y1 - k_1 * x1;
                b_2 = y1 - k_2 * x1;
                double x_1 = (2 * x0 - 2 * k_1 * b_1 + 2 * k_1 * y0) / (2 * (1 + pow(k_1, 2)));
                touch_point1 = Vec2D{x_1, k_1 * x_1 + b_1};

                double x_2 = (2 * x0 - 2 * k_2 * b_2 + 2 * k_2 * y0) / (2 * (1 + pow(k_2, 2)));
                touch_point2 = Vec2D{x_2, k_2 * x_2 + b_2};
            }

            touch_point_file << touch_point1.x << " " << touch_point1.y << "\n";
            touch_point_file << touch_point2.x << " " << touch_point2.y << "\n";

            if (j == i-1) {
                if (points[i].get_2d().lengh(touch_point1) < points[i].get_2d().lengh(touch_point2))
                    points_full.push_back(Vec3D{touch_point1, points[i].h});
                else
                    points_full.push_back(Vec3D{touch_point2, points[i].h});
            }
            else
            {
                Vec3D prev_circ_point = points_full[points_full.size() - 2];
                if ((mid_k*prev_circ_point.x + mid_b - prev_circ_point.y)*(mid_k*touch_point1.x + mid_b - touch_point1.y) < 0)
                {
                    points_full.push_back(Vec3D{touch_point1, points[i].h});
                }
                else
                    points_full.push_back(Vec3D{touch_point2, points[i].h});
            }

            if (j == i-1)
                points_full.push_back(points[i]);

        }
    }
    points_full.push_back(points[points.size() - 1]);
    points = points_full;
    center_file.close();
    touch_point_file.close();
};

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
