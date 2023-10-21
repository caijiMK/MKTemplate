#ifndef MKTemplate_Mint
#define MKTemplate_Mint

struct mint {
	static const int mod = 998244353;
	int v;

	mint() = default;
	mint(int _v): v(_v) {}
	explicit operator int() const {return v;}
	mint operator+(const mint &x) const {return v + x.v - (v + x.v < mod ? 0 : mod);}
	mint &operator+=(const mint &x) {return *this = *this + x;}
	mint operator-(const mint &x) const {return v - x.v + (v - x.v >= 0 ? 0 : mod);}
	mint &operator-=(const mint &x) {return *this = *this - x;}
	mint operator*(const mint &x) const {return (long long)v * x.v % mod;}
	mint &operator*=(const mint &x) {return *this = *this * x;}
	mint inv() const {
		mint a(*this), ans(1);
		int b(mod - 2);
		while (b) {
			if (b & 1) ans *= a;
			a *= a;
			b >>= 1;
		}
		return ans;
	}
	mint operator/(const mint &x) const {return *this * x.inv();}
	mint &operator/=(const mint &x) {return *this = *this / x;}
	mint operator-() {return mint(-v);}
};

#endif