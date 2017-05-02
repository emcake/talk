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

        fn zero() -> Self { Self::from(0u8) }
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

    fn point_within(&self, pt:&V) -> bool {

        let z = T::zero();

        (pt.x() - self.min().x()) >= z && (pt.x() - self.max().x()) <= z && (pt.y() - self.min().y()) >= (-z) && (pt.y() - self.max().y()) <= z
    }

    fn closest_bounded_point(&self, p:&V) -> V {
        let y = 
            {
                if p.y() < self.min().y() {
                    self.min().y()
                }
                else if p.y() > self.max().y() {
                    self.max().y()
                }
                else{
                    p.y()
                }
            };

        let x = 
            {
                if p.x() < self.min().x() {
                    self.min().x()
                }
                else if p.x() > self.max().x() {
                    self.max().x()
                }
                else{
                    p.x()
                }
            };

        V::new(x,y)
    }
    
}

