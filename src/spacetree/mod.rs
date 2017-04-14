use std::cmp::PartialOrd;

trait Dist<V> {
    type Output : PartialOrd;

    fn dist_between(v1 : &V, v2 : &V) -> Self::Output;
}

trait SpaceTree <'a, T : PartialOrd, V, D : Dist<V,Output=T>> {
    fn new(pts : Iterator<Item=V>) -> Self;

    fn closest_point(target : &V) -> (&V, D::Output);

    fn points_within_dist(target: &V, dist: &D::Output) -> &'a Iterator<Item=V>;
}