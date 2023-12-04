#ifndef MKTemplate_Geometry
#define MKTemplate_Geometry

#include <cmath>
#include <vector>

namespace Geometry {
	const double eps = 1e-6;
	const double pi = acos(-1);

	struct Vect {
		double x, y;

		Vect() = default;
		Vect(double _x, double _y): x(_x), y(_y) {}
		Vect operator+(const Vect &a) const {return Vect(x + a.x, y + a.y);}
		Vect operator-(const Vect &a) const {return Vect(x - a.x, y - a.y);}
		Vect operator*(const double &t) const {return Vect(x * t, y * t);}
		Vect operator/(const double &t) const {return Vect(x / t, y / t);}
		double operator*(const Vect &a) const {return x * a.x + y * a.y;}
		double operator^(const Vect &a) const {return x * a.y - y * a.x;}
		Vect &operator+=(const Vect &a) {return *this = *this + a;}
		Vect &operator-=(const Vect &a) {return *this = *this - a;}
		Vect &operator*=(const double &t) {return *this = *this * t;}
		Vect &operator/=(const double &t) {return *this = *this / t;}
		double length() {return sqrt(x * x + y * y);}
		void rotate(double rad) {
			*this = Vect(x * cos(rad) - y * sin(rad), x * sin(rad) + y * cos(rad));
			return ;
		}
	};
	using Point = Vect;
	struct Polygon: vector<Point> {
		double area() {
			double ans = 0;
			for (int i = 0; i < (int)size(); i++) ans += at(i) ^ at((i + 1) % size());
			return ans / 2;
		}
	};
	struct Line {
		Point s;
		Vect dir;

		Line() = default;
		Line(Point _s, Vect _dir): s(_s), dir(_dir) {}
	};
	struct Circle {
		Point O;
		double r;

		Circle() = default;
		Circle(Point _O, double _r): O(_O), r(_r) {}
		double sector(double rad) {return r * r * rad / 2;}
		double arc(double rad) {return (r * r * rad - r * r * sin(rad)) / 2;}
	};

	double dist(Point u, Point v) {return sqrt((u.x - v.x) * (u.x - v.x) + (u.y - v.y) * (u.y - v.y));}
	double angle(Vect a, Vect b) {return atan2(a ^ b, a * b);}
	int onSegment(Point p, Line s) {
		return fabs(s.dir ^ (p - s.s)) <= eps && (p - s.s) * (p - s.s - s.dir) <= 0;
	}
	Point intersection(Line a, Line b) {
		double t = (a.dir ^ (b.s - a.s)) / (b.dir ^ a.dir);
		return b.s + b.dir * t;
	}
	pair<Point, Point>intersection(Circle &a, Circle &b) {
		Vect u = b.O - a.O;
		double rad = acos((a.r * a.r + u.length() * u.length() - b.r * b.r) / a.r / u.length() / 2);
		Vect v = u;
		u.rotate(-rad), v.rotate(rad);
		return {u + a.O, v + a.O};
	}
	Polygon Andrew(Polygon a) {
		Polygon ans;
		sort(a.begin(), a.end(), [](const Point &u, const Point &v) {
			return u.x != v.x ? u.x < v.x : u.y < v.y;
		});
		ans.push_back(a[0]);
		for (int i = 1; i < (int)a.size(); i++) {
			while (ans.size() > 1) {
				Point tmp = ans.back();
				ans.pop_back();
				if (((tmp - ans.back()) ^ (a[i] - tmp)) > 0) {ans.push_back(tmp); break;}
			}
			ans.push_back(a[i]);
		}
		int top = ans.size();
		for (int i = (int)a.size() - 1; i >= 1; i--) {
			while ((int)ans.size() > top) {
				Point tmp = ans.back();
				ans.pop_back();
				if (((tmp - ans.back()) ^ (a[i] - tmp)) > 0) {ans.push_back(tmp); break;}
			}
			ans.push_back(a[i]);
		}
		return ans;
	}
	long long diameter(Polygon a) {
		const int n = a.size();
		if (n == 1) return 0;
		if (n == 2) return dist(a[0], a[1]);
		long long ans = 0;
		a.push_back(a[0]);
		for (int i = 0, pos = 2; i < n; i++) {
			Point u = a[i], v = a[i + 1];
			ans = max(ans, dist(u, v));
			while (((u - a[pos]) ^ (v - a[pos])) < ((u - a[pos + 1]) ^ (v - a[pos + 1])))
				pos = (pos + 1) % n;
			ans = max({ans, dist(u, a[pos]), dist(v, a[pos])});
		}
		return ans;
	}
}

#endif