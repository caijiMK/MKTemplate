#include <cassert>
#include <climits>
#include <cmath>
#include <iostream>
#include <ostream>

#include "geometry.hpp"

namespace UnitTest {
	std::ostream &operator<<(std::ostream &os, const Geometry::Vect &a) {
		return os << "(" << a.x << "," << a.y << ")";
	}

	void UnitTestSign() {
		// WARN: not work with NAN
		assert(Geometry::sign(-1) < 0);
		assert(Geometry::sign(-Geometry::eps) == 0);
		assert(Geometry::sign(0) == 0);
		assert(Geometry::sign(Geometry::eps) == 0);
		assert(Geometry::sign(1) > 0);
	}
	void UnitTestRotate() {
		assert(Geometry::rotate(Geometry::Vect(), 0) == Geometry::Vect());
		assert(Geometry::rotate(Geometry::Vect(1, 0), Geometry::pi / 2) ==
			   Geometry::Vect(0, 1));
		assert(Geometry::rotate(Geometry::Vect(1, 0), Geometry::pi) ==
			   Geometry::Vect(-1, 0));
		assert(Geometry::rotate(Geometry::Vect(1, 0), Geometry::pi * 3 / 2) ==
			   Geometry::Vect(0, -1));
		assert(Geometry::rotate(Geometry::Vect(1, 0), Geometry::pi * 2) ==
			   Geometry::Vect(1, 0));
	}
	void UnitTestDist() {
		assert(Geometry::dist(Geometry::Point(0, 0), Geometry::Point(0, 0)) ==
			   0);
		assert(Geometry::sign(std::fabs(Geometry::dist(Geometry::Point(0, 0),
													   Geometry::Point(1, 1)) -
										std::sqrt(2.0))) == 0);
	}
	void UnitTestAngle() {
		assert(Geometry::angle(Geometry::Vect(0, 0), Geometry::Vect(0, 0)) ==
			   0);
		assert(Geometry::angle(Geometry::Vect(1, 0), Geometry::Vect(1, 0)) ==
			   0);
		assert(Geometry::sign(Geometry::angle(Geometry::Vect(1, 0),
											  Geometry::Vect(0, 1)) -
							  Geometry::pi / 2) == 0);
	};
	void UnitTestProjection() {
		assert(
			Geometry::projection(Geometry::Vect(0, 0), Geometry::Vect(1, 0)) ==
			Geometry::Vect(0, 0));
		assert(
			Geometry::projection(Geometry::Vect(2, 4), Geometry::Vect(4, 0)) ==
			Geometry::Vect(2, 0));
		assert(
			Geometry::projection(Geometry::Vect(3, 0), Geometry::Vect(3, 3)) ==
			Geometry::Vect(1.5, 1.5));

		assert(Geometry::projection(Geometry::Point(0, 0),
									Geometry::Line(Geometry::Point(0, 0),
												   Geometry::Vect(1, 0))) ==
			   Geometry::Point(0, 0));
		assert(Geometry::projection(Geometry::Point(2, 4),
									Geometry::Line(Geometry::Point(-5, 0),
												   Geometry::Vect(1, 0))) ==
			   Geometry::Point(2, 0));
		assert(
			Geometry::projection(Geometry::Point(3, 0),
								 Geometry::Line(Geometry::Point(-1, -1),
												Geometry::Vect(1.5, 1.5))) ==
			Geometry::Point(1.5, 1.5));
	}
	void UnitTestReflection() {
		assert(Geometry::reflection(Geometry::Point(0, 0), Geometry::Line()) ==
			   Geometry::Point(0, 0));
		assert(Geometry::reflection(Geometry::Point(1, 0),
									Geometry::Line(Geometry::Point(0, 0),
												   Geometry::Vect(1, 0))) ==
			   Geometry::Point(-1, 0));
	}
}; // namespace UnitTest

int main() {
	UnitTest::UnitTestSign();
	UnitTest::UnitTestRotate();
	UnitTest::UnitTestDist();
	UnitTest::UnitTestAngle();
	UnitTest::UnitTestProjection();
	return 0;
}
