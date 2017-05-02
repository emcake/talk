use std::cmp::PartialOrd;
use std::fmt::Debug;

pub trait Dist<V> {
    type Output : PartialOrd + Debug;

    fn dist_between(v1 : &V, v2 : &V) -> Self::Output;
}

pub trait SpaceTree <'a, V : 'a> {
    type I : Iterator<Item=&'a V>;

    fn new<Pts : IntoIterator<Item=V>>(pts : Pts) -> Self;

    fn closest_point<'b : 'a, D : Dist<V>>(&'a self, target : &'b V) -> Option<(D::Output, &'a V)>;

    fn points_within_dist<'b : 'a, D : Dist<V>>(&'a self, target: &'b V, dist: &'b D::Output) -> Self::I;
}