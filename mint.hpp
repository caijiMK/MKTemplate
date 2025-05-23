#ifndef MKTemplate_Mint
#define MKTemplate_Mint

template<int mod>
struct modint {
	int v;

	enum {WITHMOD};
	modint(): v(0) {};
	modint(int _v): v(_v) {}
	modint(int _v, int): v((_v % mod + mod) % mod) {}
	explicit operator int() const {return v;}
	explicit operator long long() const {return v;}
	modint operator+(const modint &b) const {return v + b.v >= mod ? v + b.v - mod : v + b.v;}
	modint &operator+=(const modint &b) {v = (v + b.v >= mod ? v + b.v - mod : v + b.v); return *this;}
	modint operator-(const modint &b) const {return v - b.v < 0 ? v - b.v + mod : v - b.v;}
	modint &operator-=(const modint &b) {v = (v - b.v < 0 ? v - b.v + mod : v - b.v); return *this;}
	modint operator*(const modint &b) const {return (long long)v * b.v % mod;}
	modint &operator*=(const modint &b) {v = (long long)v * b.v % mod; return *this;}
	modint inv() const {
		int a = v, b = mod - 2, ans = 1;
		while (b) {
			if (b & 1) ans = (long long)ans * a % mod;
			a = (long long)a * a % mod;
			b >>= 1;
		}
		return ans;
	}
	modint operator/(const modint &x) const {return (long long)v * x.inv().v % mod;}
	modint &operator/=(const modint &x) {v = (long long)v * x.inv().v % mod; return *this;}
	modint operator-() {return modint(-v);}
	modint operator+() {return *this;}
	bool operator==(const modint &b) {return v == b.v;}
	bool operator!=(const modint &b) {return v != b.v;}
};
using mint = modint<998244353>;

#endif