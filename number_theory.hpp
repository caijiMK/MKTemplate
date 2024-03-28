#ifndef MKTemplate_NumberTheory
#define MKTemplate_NumberTheory

#include <vector>
using namespace std;

namespace NumberTheory {
	namespace MillerRabin {
		const int prime[] = {2, 3, 5, 7, 9, 11, 13, 17, 19, 23, 29, 31, 37};

		inline long long power(long long a, long long b, long long mod) {
			long long ans = 1;
			while (b) {
				if (b & 1) ans = (__int128)ans * a % mod;
				a = (__int128)a * a % mod;
				b >>= 1;
			}
			return ans % mod;
		}
	}
	inline int Miller_Rabin(long long n) {
		if (n == 1) return 0;
		if (n == 2) return 1;
		if (n % 2 == 0) return 0;
		long long u = n - 1, t = 0;
		while (u % 2 == 0) u /= 2, t++;
		for (int i = 0; i < 12; i++) {
			if (MillerRabin::prime[i] % n == 0) continue;
			long long x = MillerRabin::power(MillerRabin::prime[i] % n, u, n);
			if (x == 1) continue;
			int flag = 0;
			for (int j = 1; j <= t; j++) {
				if (x == n - 1) {flag = 1; break;}
				x = (__int128)x * x % n;
			}
			if (!flag) return 0;
		}
		return 1;
	}

	inline long long exgcd(long long a, long long b, long long &x, long long &y) {
		if (b == 0) {x = 1; y = 0; return a;}
		long long m = exgcd(b, a % b, y, x);
		y -= a / b * x;
		return m;
	}

	inline long long CRT(vector<pair<int, int>> vec) {
		long long ans = 0, mul = 1;
		for (auto i : vec) mul *= i.second;
		for (auto i : vec) {
			long long m = mul / i.second;
			long long x, y;
			exgcd(m, i.second, x, y);
			x = (x % i.second + i.second) % i.second;
			ans = (ans + i.first * m * x) % mul;
		}
		return ans;
	}

	inline long long exCRT(vector<pair<int, int>> vec) {
		long long ans = vec[0].first, mod = vec[0].second;
		for (int i = 1; i < (int)vec.size(); i++) {
			long long a = mod, b = vec[i].second, c = vec[i].second - ans % b;
			long long x, y;
			long long g = exgcd(a, b, x, y);
			if (c % g != 0) return -1;
			b /= g;
			x = x * (c / g) % b;
			ans += x * mod;
			mod *= b;
			ans = (ans % mod + mod) % mod;
		}
		return ans;
	}
}

#endif