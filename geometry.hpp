#ifndef MKTemplate_Geometry
#define MKTemplate_Geometry

#include <cmath>
#include <vector>
#include <algorithm>

namespace Geometry {
	using namespace std;

	const double eps = 1e-6;
	const double pi = acos(-1);

	struct Vect { // 向量 & 点
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
	struct Polygon: vector<Point> { // 多边形
		Polygon(initializer_list<Point> &&arg): vector<Point>(arg) {}
		template<typename... argT>
		Polygon(argT &&...args): vector<Point>(forward<argT>(args)...) {}

		double area() { // 多边形面积
			double ans = 0;
			for (int i = 0; i < (int)size(); i++) ans += (*this)[i] ^ (*this)[(i + 1) % size()];
			return ans / 2;
		}
	};
	struct Line { // 直线
		Point s;
		Vect dir;
		// 起点 + 向量

		Line() = default;
		Line(Point _s, Vect _dir): s(_s), dir(_dir) {}
	};
	struct Segment: Line {
		Segment() = default;
		Segment(Point _s, Vect _dir): Line(_s, _dir) {}
	};
	struct Circle { // 圆
		Point O;
		double r;
		// 圆心 + 半径

		Circle() = default;
		Circle(Point _O, double _r): O(_O), r(_r) {}
		double sector(double rad) {return r * r * rad / 2;} // 扇形面积
		double arc(double rad) {return (r * r * rad - r * r * sin(rad)) / 2;} // ?
	};

	inline int cmp0(double x) {return fabs(x) <= eps ? 0 : (x > 0 ? 1 : -1);}
	// 两点距离
	inline double dist(Point u, Point v) {return sqrt((u.x - v.x) * (u.x - v.x) + (u.y - v.y) * (u.y - v.y));}
	// 夹角
	inline double angle(Vect a, Vect b) {return atan2(a ^ b, a * b);}
	// 求向量 a 到向量 b 的投影
	inline Vect projection(Vect a, Vect b) {return b * (a * b / b.length() / b.length());}
	// 求点 a 到直线 b 的投影
	inline Point projection(Point a, Line b) {return b.s + projection(a - b.s, b.dir);}
	// 求点 a 关于直线 b 的轴对称点
	inline Point reflection(Point a, Line b) {
		Point tmp = b.s + projection(a - b.s, b.dir);
		return tmp + tmp - a;
	}
	// 判断两条直线是否平行
	inline int isParallel(Line a, Line b) {return fabs(a.dir ^ b.dir) <= eps;}
	// 判断两条直线是否垂直
	inline int isOrthogonal(Line a, Line b) {return fabs(a.dir * b.dir) <= eps;}
	// 判断点 a 是否在线段 b 上
	inline int onSegment(Point a, Segment b) {
		return fabs(b.dir ^ (a - b.s)) <= eps && cmp0((a - b.s) * (a - b.s - b.dir)) <= 0;
		// 第一条判断 a 是否在 b 所在直线上，第二条判断是否夹在两个端点中间
	}
	// 判断两条线段是否相交
	inline int isIntersecting(Segment a, Segment b) {
		return onSegment(b.s, a) || onSegment(b.s + b.dir, a) ||
			onSegment(a.s, b) || onSegment(a.s + a.dir, b) ||
			(cmp0((b.s - a.s) ^ a.dir) != cmp0((b.s + b.dir - a.s) ^ a.dir) &&
				cmp0((a.s - b.s) ^ b.dir) != cmp0((a.s + a.dir - b.s) ^ b.dir));
	}
	// 线段交点
	inline Point intersection(Segment a, Segment b) {
		double t = (a.dir ^ (b.s - a.s)) / (b.dir ^ a.dir);
		return b.s + b.dir * t;
	}
	// 点 a 到线段 b 的距离
	inline double dist(Point a, Segment b) {
		double ans = min(dist(a, b.s), dist(a, b.s + b.dir));
		Point tmp = projection(a, b);
		if (onSegment(tmp, b)) ans = min(ans, dist(a, tmp));
		return ans;
	}
	// 线段间的距离
	inline double dist(Segment a, Segment b) {
		if (isIntersecting(a, b)) return 0;
		return min({dist(a.s, b), dist(a.s + a.dir, b), dist(b.s, a), dist(b.s + b.dir, a)});
	}
	inline int isConvex(Polygon a) {
		for (int i = 0; i < (int)a.size(); i++)
			if (cmp0((a[(i + 1) % (int)a.size()] - a[i]) ^ (a[(i + 2) % (int)a.size()] - a[i])) < 0) return 0;
		return 1;
	}
	inline pair<Point, Point>intersection(Circle &a, Circle &b) { // 两个圆的交点
		Vect u = b.O - a.O;
		double rad = acos((a.r * a.r + u.length() * u.length() - b.r * b.r) / a.r / u.length() / 2);
		Vect v = u;
		u.rotate(-rad), v.rotate(rad);
		return {u + a.O, v + a.O};
	}
	inline Polygon Andrew(Polygon a) { // 凸包
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
	inline Polygon Minkowski(Polygon l, Polygon r) { // 闵可夫斯基和
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
	inline double Diameter(Polygon a) { // 旋转卡壳
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