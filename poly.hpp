#ifndef MKTemplate_Poly
#define MKTemplate_Poly

namespace Poly {
	const int mod = 998244353, G = 3, invG = 332748118;

	int power(int a, int b) {
		int ans = 1;
		while (b) {
			if (b & 1) ans = (long long)ans * a % mod;
			a = (long long)a * a % mod;
			b >>= 1;
		}
		return ans % mod;
	}

	struct poly: vector<int> {
		template<typename ... argT>
		poly(argT &&... args): vector<int>(forward<argT>(args)...) {}

		poly operator+(const poly &b) const {
			const poly &a = *this;
			poly ans(max(a.size(), b.size()));
			for (int i = 0; i < (int)ans.size(); i++)
				ans[i] = (i < (int)a.size() ? a[i] : 0) + (i < (int)b.size() ? b[i] : 0);
			return ans;
		}
		poly operator+=(const poly &b) {return *this = *this + b;}
		poly operator-(const poly &b) const {
			const poly &a = *this;
			poly ans(max(a.size(), b.size()));
			for (int i = 0; i < (int)ans.size(); i++)
				ans[i] = (i < (int)a.size() ? a[i] : 0) - (i < (int)b.size() ? b[i] : 0);
			return ans;
		}
		poly operator-=(const poly &b) {return *this = *this - b;}
		void NTT(poly &g, int flag) const {
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
						int t = (long long)w[j] * f[i + mid + j] % mod;
						f[i + mid + j] = f[i + j] - t + mod;
						f[i + j] += t;
					}
				if (mid == 1 << 10)
					for (int i = 0; i < n; i++) f[i] %= mod;
			}
			int inv = flag ? 1 : power(n, mod - 2);
			for (int i = 0; i < n; i++) g[i] = f[i] % mod * inv % mod;
			return ;
		}
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
		poly sqrt() const { // need F[0] = 1.
			poly f = *this, g;
			g.push_back(1);
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

	poly power(poly f, int b) { // need F[0] = 1.
		f = f.ln();
		for (int i = 0; i < (int)f.size(); i++) f[i] = (long long)f[i] * b % mod;
		f = f.exp();
		return f;
	}
	// poly power(poly f, int b1, int b2 = -1) {
	// 	if (b2 == -1) b2 = b1;
	// 	int n = f.size(), p = 0;
	// 	reverse(f.begin(), f.end());
	// 	while (!f.empty() && !f.back()) f.pop_back(), p++;
	// 	if (f.empty() || (long long)p * b1 >= n) return poly(n);
	// 	int v = f.back();
	// 	int inv = power(v, mod - 2);
	// 	for (int &i : f) i = (long long)i * inv % mod;
	// 	reverse(f.begin(), f.end());
	// 	f = f.ln();
	// 	for (int &i : f) i = (long long)i * b1 % mod;
	// 	f = f.exp();
	// 	reverse(f.begin(), f.end());
	// 	for (int i = 1; i <= p * b1; i++) f.push_back(0);
	// 	reverse(f.begin(), f.end());
	// 	f.resize(n);
	// 	v = power(v, b2);
	// 	for (int &i : f) i = (long long)i * v % mod;
	// 	return f;
	// }
}

#endif