#ifndef MKTemplate_Poly
#define MKTemplate_Poly

#include <chrono>
#include <random>
#include <vector>
#include <algorithm>

namespace Poly {
	using namespace std;

	const int mod = 998244353, G = 3, invG = 332748118;

	inline int power(int a, int b) {
		int ans = 1;
		while (b) {
			if (b & 1) ans = (long long)ans * a % mod;
			a = (long long)a * a % mod;
			b >>= 1;
		}
		return ans % mod;
	}

	struct poly: vector<int> {
		poly(initializer_list<int> &&arg): vector<int>(arg) {}
		template<typename... argT>
		explicit poly(argT &&...args): vector<int>(forward<argT>(args)...) {}

		poly operator+(const poly &b) const {
			const poly &a = *this;
			poly ans(max(a.size(), b.size()));
			for (int i = 0; i < (int)ans.size(); i++)
				ans[i] = ((i < (int)a.size() ? a[i] : 0) + (i < (int)b.size() ? b[i] : 0)) % mod;
			return ans;
		}
		poly operator+=(const poly &b) {return *this = *this + b;}
		poly operator-(const poly &b) const {
			const poly &a = *this;
			poly ans(max(a.size(), b.size()));
			for (int i = 0; i < (int)ans.size(); i++)
				ans[i] = ((i < (int)a.size() ? a[i] : 0) - (i < (int)b.size() ? b[i] : 0) + mod) % mod;
			return ans;
		}
		poly operator-=(const poly &b) {return *this = *this - b;}
		static void NTT(poly &g, int flag) {
			int n = g.size();
			vector<unsigned long long> f(g.begin(), g.end());
			vector<int> swp(n);
			for (int i = 0; i < n; i++) {
				swp[i] = swp[i >> 1] >> 1 | ((i & 1) * (n >> 1));
				if (i < swp[i]) std::swap(f[i], f[swp[i]]);
			}
			for (int mid = 1; mid < n; mid <<= 1) {
				int w1 = power(flag ? G : invG, (mod - 1) / mid / 2);
				vector<int> w(mid);
				w[0] = 1;
				for (int i = 1; i < mid; i++) w[i] = (long long)w[i - 1] * w1 % mod;
				for (int i = 0; i < n; i += mid << 1)
					for (int j = 0; j < mid; j++) {
						int t = f[i + mid + j] % mod * w[j] % mod;
						f[i + mid + j] = f[i + j] - t + mod;
						f[i + j] += t;
					}
				if (mid == 1 << 10)
					for (int i = 0; i < n; i++) f[i] %= mod;
			}
			int inv = flag ? 1 : power(n, mod - 2);
			for (int i = 0; i < n; i++) g[i] = f[i] % mod * inv % mod;
			return;
		}
		// 下面是基于转置原理的 NTT，相对朴素版本效率更高。
		/*
		static void NTT(poly &g, int flag) {
			int n = g.size();
			vector<int> f(g.begin(), g.end());
			if (flag) {
				for (int mid = n >> 1; mid >= 1; mid >>= 1) {
					int w1 = power(G, (mod - 1) / mid / 2);
					vector<int> w(mid);
					w[0] = 1;
					for (int i = 1; i < mid; i++) w[i] = (long long)w[i - 1] * w1 % mod;
					for (int i = 0; i < n; i += mid << 1)
						for (int j = 0; j < mid; j++) {
							int t = (long long)(f[i + j] - f[i + mid + j] + mod) * w[j] % mod;
							f[i + j] = f[i + j] + f[i + mid + j] >= mod ?
								f[i + j] + f[i + mid + j] - mod : f[i + j] + f[i + mid + j];
							f[i + mid + j] = t;
						}
				}
				for (int i = 0; i < n; i++) g[i] = f[i];
			} else {
				for (int mid = 1; mid < n; mid <<= 1) {
					int w1 = power(invG, (mod - 1) / mid / 2);
					vector<int> w(mid);
					w[0] = 1;
					for (int i = 1; i < mid; i++) w[i] = (long long)w[i - 1] * w1 % mod;
					for (int i = 0; i < n; i += mid << 1)
						for (int j = 0; j < mid; j++) {
							int t = (long long)w[j] * f[i + mid + j] % mod;
							f[i + mid + j] = f[i + j] - t < 0 ? f[i + j] - t + mod : f[i + j] - t;
							f[i + j] = f[i + j] + t >= mod ? f[i + j] + t - mod : f[i + j] + t;
						}
				}
				int inv = power(n, mod - 2);
				for (int i = 0; i < n; i++) g[i] = (long long)f[i] * inv % mod;
			}
			return;
		}
		*/
		poly operator*(poly b) const {
			poly a(*this);
			int n = 1, len = (int)(a.size() + b.size()) - 1;
			while (n < len) n <<= 1;
			a.resize(n), b.resize(n);
			NTT(a, 1), NTT(b, 1);
			poly c(n);
			for (int i = 0; i < n; i++) c[i] = (long long)a[i] * b[i] % mod;
			NTT(c, 0);
			c.resize(len);
			return c;
		}
		poly operator*=(const poly &b) {return *this = *this * b;}
		poly inv() const {
			poly f = *this, g;
			g.push_back(power(f[0], mod - 2));
			int n = 1;
			while (n < (int)f.size()) n <<= 1;
			f.resize(n << 1);
			for (int len = 2; len <= n; len <<= 1) {
				poly tmp(len), ff(len << 1);
				for (int i = 0; i < len >> 1; i++) tmp[i] = g[i] * 2 % mod;
				for (int i = 0; i < len; i++) ff[i] = f[i];
				g.resize(len << 1);
				NTT(g, 1), NTT(ff, 1);
				for (int i = 0; i < len << 1; i++) g[i] = (long long)g[i] * g[i] % mod * ff[i] % mod;
				NTT(g, 0);
				g.resize(len);
				for (int i = 0; i < len; i++) g[i] = (tmp[i] - g[i] + mod) % mod;
			}
			g.resize(size());
			return g;
		}
		poly sqrt() const { // need F[0] != 0.
			static auto Cipolla = [](int n) -> int {
				if (n == 1) return 1;
				static auto Legendre = [](int n) -> int {return power(n, (mod - 1) / 2);};
				// assert(Legendre(n) == 1);
				static mt19937 gen(chrono::system_clock::now().time_since_epoch().count());
				int r;
				do r = gen() % mod;
				while (Legendre(((long long)r * r - n + mod) % mod) != mod - 1);
				static int sqw;
				sqw = ((long long)r * r - n + mod) % mod;
				struct complex {
					int r, i;

					complex() = default;
					complex(int _r, int _i): r(_r), i(_i) {}
					complex operator*(const complex &b) const {
						return complex(((long long)r * b.r + (long long)i * b.i % mod * sqw) % mod,
							((long long)r * b.i + (long long)i * b.r) % mod);
					}
				};
				static auto power = [](complex a, int b) -> complex {
					complex ans(1, 0);
					while (b) {
						if (b & 1) ans = ans * a;
						a = a * a;
						b >>= 1;
					}
					return ans;
				};
				int ans = power(complex(r, 1), (mod + 1) / 2).r;
				ans = min(ans, mod - ans);
				return ans;
			};
			poly f = *this, g;
			g.push_back(Cipolla(f[0]));
			int n = 1;
			while (n < (int)f.size()) n <<= 1;
			f.resize(n << 1);
			for (int len = 2; len <= n; len <<= 1) {
				poly tmp(len), ff(len << 1);
				for (int i = 0; i < len >> 1; i++) tmp[i] = g[i] * 2 % mod;
				for (int i = 0; i < len; i++) ff[i] = f[i];
				g.resize(len << 1);
				NTT(g, 1);
				for (int i = 0; i < len << 1; i++) g[i] = (long long)g[i] * g[i] % mod;
				NTT(g, 0);
				g += ff;
				g *= tmp.inv();
				g.resize(len);
			}
			g.resize(size());
			return g;
		}
		poly derivative() const {
			poly f(*this);
			for (int i = 1; i < (int)f.size(); i++) f[i - 1] = (long long)f[i] * i % mod;
			f.pop_back();
			return f;
		}
		poly integral() const {
			poly f(*this);
			f.push_back(0);
			for (int i = f.size() - 1; i >= 1; i--) f[i] = (long long)f[i - 1] * power(i, mod - 2) % mod;
			f[0] = 0;
			return f;
		}
		poly ln() const {
			poly f((derivative() * inv()).integral());
			f.resize(size()); 
			return f;
		}
		poly exp() const { // need F[0] = 0.
			poly f(*this), g;
			g.push_back(1);
			int n = 1;
			while (n < (int)size()) n <<= 1;
			f.resize(n);
			for (int len = 2; len <= n; len <<= 1) {
				poly tmp(g);
				g.resize(len);
				g = g.ln();
				for (int i = 0; i < len; i++) g[i] = (f[i] - g[i] + mod) % mod;
				g[0] = (g[0] + 1) % mod;
				g *= tmp;
				g.resize(len);
			}
			g.resize(size());
			return g;
		}
	};

	inline poly power(poly f, int b) { // need F[0] = 1.
		f = f.ln();
		for (int i = 0; i < (int)f.size(); i++) f[i] = (long long)f[i] * b % mod;
		f = f.exp();
		return f;
	}
	/*
	不要求 F[0] = 1
	b1 为指数 mod m，b2 为指数 mod (phi(m) = mod - 1)
	如果 b1 传入之前被取模了记得特判 a[0] = 0 且指数 > n 的情况，此时多项式每一位系数都是 0
	*/
	/*
	inline poly power(poly f, int b1, int b2 = -1) {
		if (b2 == -1) b2 = b1;
		int n = f.size(), p = 0;
		reverse(f.begin(), f.end());
		while (!f.empty() && !f.back()) f.pop_back(), p++;
		if (f.empty() || (long long)p * b1 >= n) return poly(n);
		int v = f.back();
		int inv = power(v, mod - 2);
		for (int &i : f) i = (long long)i * inv % mod;
		reverse(f.begin(), f.end());
		f = f.ln();
		for (int &i : f) i = (long long)i * b1 % mod;
		f = f.exp();
		reverse(f.begin(), f.end());
		for (int i = 1; i <= p * b1; i++) f.push_back(0);
		reverse(f.begin(), f.end());
		f.resize(n);
		v = power(v, b2);
		for (int &i : f) i = (long long)i * v % mod;
		return f;
	}
	*/
	inline pair<poly, poly> div(poly f, poly g) {
		int n = f.size() - 1, m = g.size() - 1;
		reverse(f.begin(), f.end()), reverse(g.begin(), g.end());
		g.resize(n - m + 1);
		poly ans = f * g.inv();
		ans.resize(n - m + 1);
		reverse(ans.begin(), ans.end());
		return {ans, f - g * ans};
	}
	inline vector<int> eval(poly f, vector<int> vec) {
		static vector<poly> mul;
		int n = max(f.size(), vec.size()), m = vec.size();
		mul.resize(4 * n + 5);
		f.resize(n), vec.resize(n);
		auto init = [&vec](auto &&self, int l, int r, int now = 1) -> void {
			if (l == r) {mul[now] = poly({1, (mod - vec[l]) % mod}); return;}
			int mid = (l + r) / 2, ls = now * 2, rs = now * 2 + 1;
			self(self, l, mid, ls), self(self, mid + 1, r, rs);
			mul[now] = mul[ls] * mul[rs];
			return;
		};
		init(init, 0, n - 1);
		auto mult = [](poly f, poly g) -> poly {
			int n = f.size(), m = g.size();
			reverse(g.begin(), g.end());
			g *= f;
			for (int i = 0; i < n; i++) f[i] = g[i + m - 1];
			return f;
		};
		auto solve = [&vec, &mult](auto &&self, int l, int r, poly f, int now = 1) -> void {
			f.resize(r - l + 1);
			if (l == r) {vec[l] = f[0]; return;}
			int mid = (l + r) / 2, ls = now * 2, rs = now * 2 + 1;
			self(self, l, mid, mult(f, mul[rs]), ls), self(self, mid + 1, r, mult(f, mul[ls]), rs);
			return;
		};
		solve(solve, 0, n - 1, mult(f, mul[1].inv()));
		vec.resize(m);
		return vec;
	}
	inline poly interpolation(vector<pair<int, int>> vec) {
		static vector<poly> mul;
		int n = vec.size();
		mul.resize(4 * n + 5);
		poly x(n), y(n);
		for (int i = 0; i < n; i++) x[i] = vec[i].first, y[i] = vec[i].second;
		auto init = [&x](auto &&self, int l, int r, int now = 1) -> void {
			if (l == r) {mul[now] = poly({(mod - x[l]) % mod, 1}); return;}
			int mid = (l + r) / 2, ls = now * 2, rs = now * 2 + 1;
			self(self, l, mid, ls), self(self, mid + 1, r, rs);
			mul[now] = mul[ls] * mul[rs];
			return;
		};
		init(init, 0, n - 1);
		vector<int> val = eval(mul[1].derivative(), x);
		auto solve = [&y, &val](auto &&self, int l, int r, int now = 1) -> poly {
			if (l == r) return poly({int((long long)y[l] * power(val[l], mod - 2) % mod)});
			int mid = (l + r) / 2, ls = now * 2, rs = now * 2 + 1;
			return self(self, l, mid, ls) * mul[rs] + self(self, mid + 1, r, rs) * mul[ls];
		};
		return solve(solve, 0, n - 1);
	}
}

#endif