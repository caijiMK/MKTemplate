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
		long long operator*(const Vect &a) const {return (long long)x * a.x + (long long)y * a.y;}
		long long operator^(const Vect &a) const {return (long long)x * a.y - (long long)y * a.x;}
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
	using polygon = vector<Vect>;
	struct Line {
		Point s;
		Vect dir;

		Line() = default;
		Line(Point _s, Vect _dir): s(_s), dir(_dir) {}
	};

	long long dist(Point u, Point v) {
		return (long long)(u.x - v.x) * (u.x - v.x) + (long long)(u.y - v.y) * (u.y - v.y);
	}
	int onSegment(Point p, Line s) {
		return fabs(s.dir ^ (p - s.s)) <= eps && (p - s.s) * (p - s.s - s.dir) <= 0;
	}
	Point intersection(Line a, Line b) {
		Vect p = b.s - a.s;
		double t = (a.dir ^ p) / (b.dir ^ a.dir);
		return b.s + b.dir * t;
	}
	double angle(Vect a, Vect b) {return acos(min(max(a * b / a.length() / b.length(), -1.0), 1.0));}
	polygon Andrew(polygon a) {
		polygon ans;
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
	long long diameter(polygon a) {
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