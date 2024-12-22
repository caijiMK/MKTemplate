#ifndef MKTemplate_Geometry
#define MKTemplate_Geometry

#include <cmath>
#include <queue>
#include <chrono>
#include <random>
#include <vector>
#include <algorithm>

namespace Geometry {
	using namespace std;

	const double eps = 1e-9;
	const double pi = acos(-1);

	// 判断 x < 0 或 x = 0 或 x > 0
	inline int sign(double x) {return fabs(x) <= eps ? 0 : x > 0 ? 1 : -1;}

	// 向量 & 点
	struct Vect {
		double x, y;

		Vect() = default;
		Vect(double _x, double _y): x(_x), y(_y) {}
		Vect operator+() const {return *this;}
		Vect operator-() const {return Vect(-x, -y);}
		Vect operator+(const Vect &b) const {return Vect(x + b.x, y + b.y);}
		Vect operator-(const Vect &b) const {return Vect(x - b.x, y - b.y);}
		Vect operator*(const double &t) const {return Vect(x * t, y * t);}
		Vect operator/(const double &t) const {return Vect(x / t, y / t);}
		double operator*(const Vect &b) const {return x * b.x + y * b.y;}
		double operator^(const Vect &b) const {return x * b.y - y * b.x;}
		Vect &operator+=(const Vect &b) {return *this = *this + b;}
		Vect &operator-=(const Vect &b) {return *this = *this - b;}
		Vect &operator*=(const double &t) {return *this = *this * t;}
		Vect &operator/=(const double &t) {return *this = *this / t;}
		bool operator==(const Vect &b) {return !sign(x - b.x) && !sign(y - b.y);}
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
		double area() {return r * r * pi;}
	};

	// 两点的中点
	inline Point middle(Point a, Point b) {return Point((a.x + b.x) / 2, (a.y + b.y) / 2);}
	// 向量 u 逆时针旋转 rad 后得到的向量
	inline Vect rotate(Vect a, double rad) {
		return Vect(a.x * cos(rad) - a.y * sin(rad), a.x * sin(rad) + a.y * cos(rad));
	}
	// 两点距离
	inline double dist(Point a, Point b) {return hypot(a.x - b.x, a.y - b.y);}
	// 判断点是否在直线左侧
	inline int isLeft(Point a, Line b) {return sign(b.dir ^ (a - b.s)) > 0;}
	// 判断点是否在直线右侧
	inline int isRight(Point a, Line b) {return sign(b.dir ^ (a - b.s)) < 0;}
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
	// 判断两条直线是否相等
	inline int isSame(Line a, Line b) {return isParallel(a, b) && !sign(b.dir ^ (a.s - b.s));}
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
	inline Polygon andrew(Polygon a) {
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
	inline Polygon minkowski(Polygon a, Polygon b) {
		Polygon diffa, diffb;
		int siza = a.size(), sizb = b.size();
		for (int i = 1; i < siza; i++) diffa.push_back(a[i] - a[i - 1]);
		for (int i = 1; i < sizb; i++) diffb.push_back(b[i] - b[i - 1]);
		Polygon ans;
		ans.push_back(a[0] + b[0]);
		int pl = 0, pr = 0;
		while (pl < siza - 1 && pr < sizb - 1)
			if (sign(diffa[pl] ^ diffb[pr]) < 0) ans.push_back(ans.back() + diffa[pl++]);
			else ans.push_back(ans.back() + diffb[pr++]);
		while (pl < siza - 1) ans.push_back(ans.back() + diffa[pl++]);
		while (pr < sizb - 1) ans.push_back(ans.back() + diffb[pr++]);
		return ans;
	}
	// 旋转卡壳求凸包直径
	inline double diameter(Polygon a) {
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
	// 半平面交
	inline pair<Polygon, vector<Line>> halfplane(vector<Line> a) {
		deque<Line> ln;
		deque<Point> pt;
		sort(a.begin(), a.end(), [](const Line &x, const Line &y) -> int {
			double radx = atan2(x.dir.y, x.dir.x), rady = atan2(y.dir.y, y.dir.x);
			if (sign(radx - rady)) return sign(radx - rady) == -1;
			else return isLeft(x.s, y);
		});
		for (int i = 0; i < (int)a.size(); i++)
			if (!i || !isParallel(a[i - 1], a[i])) {
				while (!pt.empty() && isRight(pt.back(), a[i])) pt.pop_back(), ln.pop_back();
				while (!pt.empty() && isRight(pt.front(), a[i])) pt.pop_front(), ln.pop_front();
				if (ln.size() == 1 &&
					sign(atan2(a[i].dir.y, a[i].dir.x) - atan2(ln.back().dir.y, ln.back().dir.x) - pi) >= 0)
					return {Polygon(), vector<Line>()};
				if (!ln.empty()) {
					if (!sign(atan2(a[i].dir.y, a[i].dir.x) - atan2(ln.back().dir.y, ln.back().dir.x) - pi))
						return {Polygon(), vector<Line>()};
					pt.push_back(intersection(ln.back(), a[i]));
				}
				ln.push_back(a[i]);
			}
		while (pt.size() >= 2 && isRight(pt.back(), ln.front())) pt.pop_back(), ln.pop_back();
		if (ln.size() <= 2) return {Polygon(), vector<Line>()};
		pt.push_back(intersection(ln.front(), ln.back()));
		return {Polygon(pt.begin(), pt.end()), vector<Line>(ln.begin(), ln.end())};
	}
	// 扇形面积
	inline double sector(double r, double rad) {return r * r * rad / 2;}
	// 弧度为 rad 的弓形的面积
	inline double arch(double r, double rad) {return sector(r, rad) - r * r * sin(rad) / 2;}
	// 三角形内切圆
	inline Circle incircle(Point a, Point b, Point c) {
		Line u(a, rotate(b - a, angle(b - a, c - a) / 2)), v(b, rotate(a - b, angle(a - b, c - b) / 2));
		Circle ans;
		ans.O = intersection(u, v);
		ans.r = dist(ans.O, Line(a, b - a));
		return ans;
	}
	// 三角形外接圆
	inline Circle circumcircle(Point a, Point b, Point c) {
		Line u(middle(a, b), rotate(b - a, pi / 2)), v(middle(a, c), rotate(c - a, pi / 2));
		Point O = intersection(u, v);
		return Circle(O, dist(O, a));
	}
	// 最小圆覆盖
	inline Circle minimumCircle(Polygon a) {
		mt19937 gen(chrono::system_clock::now().time_since_epoch().count());
		int n = a.size();
		shuffle(a.begin(), a.end(), gen);
		Circle ans(Point(0, 0), 0);
		for (int i = 0; i < n; i++)
			if (sign(dist(ans.O, a[i]) - ans.r) > 0) {
				ans = Circle(a[i], 0);
				for (int j = 0; j < i; j++)
					if (sign(dist(ans.O, a[j]) - ans.r) > 0) {
						ans = Circle(middle(a[i], a[j]), dist(a[i], a[j]) / 2);
						for (int k = 0; k < j; k++)
							if (sign(dist(ans.O, a[k]) - ans.r) > 0)
								ans = circumcircle(a[i], a[j], a[k]);
					}
			}
		return ans;
	}
	// 圆与直线的交点
	inline pair<Point, Point> intersection(Circle a, Line b) {
		Point p = projection(a.O, b);
		double dis = dist(a.O, p);
		dis = sqrt(a.r * a.r - dis * dis);
		return {p - b.dir.normal() * dis, p + b.dir.normal() * dis};
	}
	// 两个圆的交点
	inline pair<Point, Point> intersection(Circle a, Circle b) {
		Vect u = b.O - a.O;
		double rad = acos((a.r * a.r + u.length() * u.length() - b.r * b.r) / a.r / u.length() / 2);
		u = u.normal() * a.r;
		return {a.O + rotate(u, rad), a.O + rotate(u, -rad)};
	}
	// 点到圆的切点
	inline pair<Point, Point> tangent(Point a, Circle b) {
		double rad = asin(b.r / dist(b.O, a));
		double dis = dist(b.O, a);
		dis = sqrt(dis * dis - b.r * b.r);
		Vect u = (b.O - a).normal() * dis;
		return {a + rotate(u, -rad), a + rotate(u, rad)};
	}
	// 两个圆的公切线，返回切点
	inline vector<pair<Point, Point>> commonTangent(Circle a, Circle b) {
		vector<pair<Point, Point>> ans;
		double dis = dist(a.O, b.O);
		if (sign(dis - fabs(a.r - b.r)) < 0) return ans;
		else if (sign(dis - fabs(a.r - b.r)) == 0) {
			ans.push_back(intersection(a, b));
			return ans;
		}
		int flag = 0;
		if (a.r < b.r) swap(a, b), flag = 1;
		Vect u = (b.O - a.O).normal();
		double rad = acos((a.r - b.r) / dis);
		ans.emplace_back(a.O + rotate(u, rad) * a.r, b.O + rotate(-u, rad - pi) * b.r);
		ans.emplace_back(a.O + rotate(u, -rad) * a.r, b.O + rotate(-u, pi - rad) * b.r);
		if (sign(dis - a.r - b.r) >= 0) {
			rad = acos((a.r + b.r) / dis);
			if (sign(rad)) {
				ans.emplace_back(a.O + rotate(u, rad) * a.r, b.O + rotate(-u, rad) * b.r);
				ans.emplace_back(a.O + rotate(u, -rad) * a.r, b.O + rotate(-u, -rad) * b.r);
			} else ans.emplace_back(a.O + u * a.r, b.O - u * b.r);
		}
		if (flag)
			for (auto &i : ans) swap(i.first, i.second);
		return ans;
	}
	// 两个圆的交的面积
	inline double intersectionArea(Circle a, Circle b) {
		int state = commonTangent(a, b).size();
		if (state <= 1) return min(a.area(), b.area());
		if (state >= 3) return 0;
		pair<Point, Point> tmp = intersection(a, b);
		double rada = fmod(angle(tmp.second - a.O, tmp.first - a.O) + 2 * pi, 2 * pi);
		double radb = fmod(angle(tmp.first - b.O, tmp.second - b.O) + 2 * pi, 2 * pi);
		return arch(a.r, rada) + arch(b.r, radb);
	}
}

#endif