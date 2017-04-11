use std::ops::Sub;
use std::ops::Mul;
use std::ops::Add;
use std::ops::Neg;
use std::ops::Div;
use std::cmp::PartialOrd;

use std::fmt::Debug;

pub trait Vectorable : Copy + Debug
        + Sub<Output=Self>
        + Mul<Output=Self>
        + Add<Output=Self> 
        + Neg<Output=Self>
        + Div<Output=Self>
        + From<u8>
        + PartialOrd
    {
        fn sqrt(&self) -> Self;

        fn max_value() -> Self;
        fn min_value() -> Self;

        fn two() -> Self { Self::from(2u8) }
    }

pub trait Vector2<T> : Copy + Debug where T : Vectorable {
    fn x(&self) -> T;
    fn y(&self) -> T;

    fn new(x:T, y:T) -> Self;

    fn dist_sq(&self, other:&Self) -> T {
        let dx : T = other.x() - self.x();
        let dy : T = other.y() - self.y();

        dx*dx + dy*dy
    }
}

use std::marker::Sized;

pub trait Bound2<V,T> : Debug + Copy
    where V : Vector2<T>, T : Vectorable, Self : Sized
{
    fn min(&self) -> &V;
    fn max(&self) -> &V;

    fn new(min:V, max:V) -> Self;

    fn split(&self, x:T, y:T) -> [Self;4] {
        [
            Self::new(V::new(self.min().x(), self.min().y()), V::new(x, y)),
            Self::new(V::new(x, self.min().y()), V::new(self.max().x(), y)),
            Self::new(V::new(self.min().x(), y), V::new(x, self.max().y())),
            Self::new(V::new(x, y), V::new(self.max().x(), self.max().y())),
        ]
    }

    fn point_within(&self, pt:&V, d:T) -> bool {

        (pt.x() - self.min().x()) >= (-d) && (pt.x() - self.max().x()) <= d && (pt.y() - self.min().y()) >= (-d) && (pt.y() - self.max().y()) <= d
    }
}

