#ifndef MKTemplate_Geometry
#define MKTemplate_Geometry

#include <cmath>
#include <vector>
#include <algorithm>
using namespace std;

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

	inline double dist(Point u, Point v) {return sqrt((u.x - v.x) * (u.x - v.x) + (u.y - v.y) * (u.y - v.y));}
	inline double angle(Vect a, Vect b) {return atan2(a ^ b, a * b);}
	inline int onSegment(Point p, Line s) {
		return fabs(s.dir ^ (p - s.s)) <= eps && (p - s.s) * (p - s.s - s.dir) <= 0;
	}
	inline Point intersection(Line a, Line b) {
		double t = (a.dir ^ (b.s - a.s)) / (b.dir ^ a.dir);
		return b.s + b.dir * t;
	}
	inline pair<Point, Point>intersection(Circle &a, Circle &b) {
		Vect u = b.O - a.O;
		double rad = acos((a.r * a.r + u.length() * u.length() - b.r * b.r) / a.r / u.length() / 2);
		Vect v = u;
		u.rotate(-rad), v.rotate(rad);
		return {u + a.O, v + a.O};
	}
	inline Polygon Andrew(Polygon a) {
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
	inline Polygon Minkowski(Polygon l, Polygon r) {
		Polygon diffl, diffr;
		int sizl = l.size(), sizr = r.size();
		for (int i = 1; i < sizl; i++) diffl.push_back(l[i] - l[i - 1]);
		for (int i = 1; i < sizr; i++) diffr.push_back(r[i] - r[i - 1]);
		Polygon ans;
		ans.push_back(l[0] + r[0]);
		int pl = 0, pr = 0;
		while (pl < sizl - 1 && pr < sizr - 1)
			if ((diffl[pl] ^ diffr[pr]) < 0) ans.push_back(ans.back() + diffl[pl++]);
			else ans.push_back(ans.back() + diffr[pr++]);
		while (pl < sizl - 1) ans.push_back(ans.back() + diffl[pl++]);
		while (pr < sizr - 1) ans.push_back(ans.back() + diffr[pr++]);
		return ans;
	}
	inline double Diameter(Polygon a) {
		const int n = a.size();
		if (n == 1) return 0;
		if (n == 2) return dist(a[0], a[1]);
		double ans = 0;
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