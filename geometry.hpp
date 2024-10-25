#ifndef MKTemplate_Geometry
#define MKTemplate_Geometry

#include <cmath>
#include <vector>
#include <algorithm>

namespace Geometry {
	using namespace std;

	const double eps = 1e-6;
	const double pi = acos(-1);

	// 向量 & 点
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
		double length() const {return sqrt(x * x + y * y);}
		Vect normal() const {return *this / length();}
	};
	using Point = Vect;
	// 多边形
	struct Polygon: vector<Point> {
		Polygon(initializer_list<Point> &&arg): vector<Point>(arg) {}
		template<typename... argT>
		explicit Polygon(argT &&...args): vector<Point>(forward<argT>(args)...) {}

		// 多边形面积
		double area() const {
			double ans = 0;
			for (int i = 0; i < (int)size(); i++) ans += (*this)[i] ^ (*this)[(i + 1) % size()];
			return ans / 2;
		}
	};
	// 直线
	struct Line {
		// 表示方式为起点 + 向量
		Point s;
		Vect dir;

		Line() = default;
		Line(Point _s, Vect _dir): s(_s), dir(_dir) {}
	};
	// 线段
	struct Segment: Line {
		Segment() = default;
		Segment(Point _s, Vect _dir): Line(_s, _dir) {}
	};
	// 圆
	struct Circle {
		// 表示方式为圆心 + 半径
		Point O;
		double r;

		Circle() = default;
		Circle(Point _O, double _r): O(_O), r(_r) {}
	};

	// 判断 x < 0 或 x = 0 或 x > 0
	inline int sign(double x) {return fabs(x) <= eps ? 0 : (x > 0 ? 1 : -1);}
	// 向量 u 逆时针旋转 rad 后得到的向量
	inline Vect rotate(Vect u, double rad) {
		return Vect(u.x * cos(rad) - u.y * sin(rad), u.x * sin(rad) + u.y * cos(rad));
	}
	// 两点距离
	inline double dist(Point u, Point v) {return sqrt((u.x - v.x) * (u.x - v.x) + (u.y - v.y) * (u.y - v.y));}
	// 两个向量的夹角
	inline double angle(Vect a, Vect b) {return atan2(a ^ b, a * b);}
	// 向量 a 到向量 b 的投影
	inline Vect projection(Vect a, Vect b) {return b * (a * b / b.length() / b.length());}
	// 点 a 到直线 b 的投影
	inline Point projection(Point a, Line b) {return b.s + projection(a - b.s, b.dir);}
	// 点 a 关于直线 b 的轴对称点
	inline Point reflection(Point a, Line b) {return b.s + projection(a - b.s, b.dir) * 2 - a;}
	// 判断两条直线是否平行
	inline int isParallel(Line a, Line b) {return !sign(a.dir ^ b.dir);}
	// 判断两条直线是否相交
	inline int isIntersecting(Line a, Line b) {return !isParallel(a, b);}
	// 判断两条直线是否垂直
	inline int isOrthogonal(Line a, Line b) {return !sign(a.dir * b.dir);}
	// 判断点 a 是否在直线 b 上
	inline int isOnLine(Point a, Line b) {return !sign(b.dir ^ (a - b.s));}
	// 判断点 a 是否在线段 b 上
	inline int isOnSegment(Point a, Segment b) {
		return isOnLine(a, b) && sign((a - b.s) * (a - b.s - b.dir)) <= 0;
	}
	// 判断直线 a 与线段 b 是否相交
	inline int isIntersecting(Line a, Segment b) {
		int va = sign((b.s - a.s) ^ a.dir), vb = sign((b.s + b.dir - a.s) ^ a.dir);
		return !va || !vb || va != vb;
	}
	// 判断两条线段是否相交
	inline int isIntersecting(Segment a, Segment b) {
		return isOnSegment(b.s, a) || isOnSegment(b.s + b.dir, a) ||
			isOnSegment(a.s, b) || isOnSegment(a.s + a.dir, b) ||
			(sign((b.s - a.s) ^ a.dir) != sign((b.s + b.dir - a.s) ^ a.dir) &&
				sign((a.s - b.s) ^ b.dir) != sign((a.s + a.dir - b.s) ^ b.dir));
	}
	// 直线交点
	inline Point intersection(Line a, Line b) {return b.s + b.dir * (a.dir ^ (b.s - a.s)) / (b.dir ^ a.dir);}
	// 点 a 到直线 b 的距离
	inline double dist(Point a, Line b) {return dist(a, projection(a, b));}
	// 点 a 到线段 b 的距离
	inline double dist(Point a, Segment b) {
		double ans = min(dist(a, b.s), dist(a, b.s + b.dir));
		Point tmp = projection(a, b);
		if (isOnSegment(tmp, b)) ans = min(ans, dist(a, tmp));
		return ans;
	}
	// 线段间的距离
	inline double dist(Segment a, Segment b) {
		if (isIntersecting(a, b)) return 0;
		return min({dist(a.s, b), dist(a.s + a.dir, b), dist(b.s, a), dist(b.s + b.dir, a)});
	}
	// 判断多边形是否是凸包
	inline int isConvex(const Polygon &a) {
		int n = a.size();
		for (int i = 0; i < n; i++)
			if (sign((a[(i + 1) % n] - a[i]) ^ (a[(i + 2) % n] - a[i])) < 0) return 0;
		return 1;
	}
	// 判断点 a 与多边形 b 的位置关系，返回值：0 - 在多边形外, 1 - 在多边形边上, 2 - 在多边形内部
	inline int Containment(Point a, const Polygon &b) {
		int n = b.size(), val = 0;
		for (int i = 0; i < n; i++) {
			Segment s = Segment(b[i], b[(i + 1) % n] - b[i]);
			if (isOnSegment(a, s)) return 1;
			int t = sign((b[(i + 1) % n] - b[i]) ^ (a - b[i]));
			int va = sign(b[i].y - a.y), vb = sign(b[(i + 1) % n].y - a.y);
			if (t > 0 && va <= 0 && vb > 0) val++;
			if (t < 0 && vb <= 0 && va > 0) val--;
		}
		return val ? 2 : 0;
	}
	// 点集的凸包
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
				if (sign((tmp - ans.back()) ^ (a[i] - tmp)) > 0) {ans.push_back(tmp); break;}
			}
			ans.push_back(a[i]);
		}
		int top = ans.size();
		for (int i = (int)a.size() - 2; i >= 0; i--) {
			while ((int)ans.size() > top) {
				Point tmp = ans.back();
				ans.pop_back();
				if (sign((tmp - ans.back()) ^ (a[i] - tmp)) > 0) {ans.push_back(tmp); break;}
			}
			ans.push_back(a[i]);
		}
		ans.pop_back();
		return ans;
	}
	// 两个凸包的闵可夫斯基和
	inline Polygon Minkowski(const Polygon &l, const Polygon &r) {
		Polygon diffl, diffr;
		int sizl = l.size(), sizr = r.size();
		for (int i = 1; i < sizl; i++) diffl.push_back(l[i] - l[i - 1]);
		for (int i = 1; i < sizr; i++) diffr.push_back(r[i] - r[i - 1]);
		Polygon ans;
		ans.push_back(l[0] + r[0]);
		int pl = 0, pr = 0;
		while (pl < sizl - 1 && pr < sizr - 1)
			if (sign(diffl[pl] ^ diffr[pr]) < 0) ans.push_back(ans.back() + diffl[pl++]);
			else ans.push_back(ans.back() + diffr[pr++]);
		while (pl < sizl - 1) ans.push_back(ans.back() + diffl[pl++]);
		while (pr < sizr - 1) ans.push_back(ans.back() + diffr[pr++]);
		return ans;
	}
	// 旋转卡壳求凸包直径
	inline double Diameter(Polygon a) {
		const int n = a.size();
		if (n == 1) return 0;
		if (n == 2) return dist(a[0], a[1]);
		double ans = 0;
		a.push_back(a[0]);
		for (int i = 0, pos = 2; i < n; i++) {
			Point u = a[i], v = a[i + 1];
			ans = max(ans, dist(u, v));
			while (sign(((u - a[pos + 1]) ^ (v - a[pos + 1])) - ((u - a[pos]) ^ (v - a[pos]))) >= 0)
				pos = (pos + 1) % n;
			ans = max({ans, dist(u, a[pos]), dist(v, a[pos])});
		}
		return ans;
	}
	// 扇形面积
	inline double sector(double r, double rad) {return r * r * rad / 2;}
	// 弧度为 rad 的弓形的面积
	inline double arch(double r, double rad) {return sector(r, rad) - r * r * sin(rad) / 2;}
	// 两个圆的公切线个数（判断两个圆的位置关系）
	inline int countOfCommonTangentLines(Circle a, Circle b) {
		if (sign(dist(a.O, b.O) - a.r - b.r) > 0) return 4;
		else if (sign(dist(a.O, b.O) - a.r - b.r) == 0) return 3;
		else if (sign(dist(a.O, b.O) - fabs(a.r - b.r)) == 0) return 1;
		else if (sign(dist(a.O, b.O) - fabs(a.r - b.r)) < 0) return 0;
		else return 2;
	}
	// 两个圆的交点
	inline pair<Point, Point>intersection(Circle a, Circle b) {
		Vect u = b.O - a.O;
		double rad = acos((a.r * a.r + u.length() * u.length() - b.r * b.r) / a.r / u.length() / 2);
		return {rotate(u, -rad) + a.O, rotate(u, rad) + a.O};
	}
}

#endif