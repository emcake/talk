use vec::Vector2;
use vec::Bound2;
use vec::Vectorable;

use std::marker::PhantomData;

use spacetree::SpaceTree;
use spacetree::Dist;

use std::fmt::Debug;

pub trait BSPBoundSplit<V> : Debug {
    fn pick_search_order(&self, v : &V) -> [usize;4];
}

pub trait BSPBound<V> : Sized + Debug + Copy {
    type S : BSPBoundSplit<V>;
    fn new<P : Iterator<Item=V>>(pts:P) -> Self;

    fn clamp(&self, v : &V) -> V;

    fn inside(&self, v : &V) -> bool;

    fn carve(&self) -> (Self::S,[Self;4]);
}

#[derive(Debug)]
pub struct Quadtree<B : BSPBound<V> ,V : Vector2<T>,T : Vectorable> (Box<QuadtreeInner<B,V,T>>,B);

#[derive(Debug)]
pub enum QuadtreeInner<B,V,T> where B : BSPBound<V>, V : Vector2<T>, T : Vectorable {
    Branch {
        split : B::S,
        children : [Quadtree<B,V,T>; 4]
    },
    Leaf {
        points : Vec<V>,
        phantom : PhantomData<T>
    }
}

impl<B,V,T> Quadtree<B,V,T> where B : BSPBound<V>, V : Vector2<T>, T : Vectorable 
{
    fn make_inner<P : IntoIterator<Item=V>>(pts : P, bound : B, itr : i32) -> Self {
        let mut ps : Vec<V> = pts.into_iter().collect();
        if ps.len() <= 1 || itr >= 10 {
            Quadtree(Box::new(QuadtreeInner::Leaf {points : ps, phantom : PhantomData}), bound)
        }
        else {

            let (split, bounds) = bound.carve();

            use arrayvec::ArrayVec;

            let kids : [Self;4] = 
                {
                    let output : ArrayVec<[Self;4]> = 
                        bounds.into_iter().scan(0, |st, b| {
                            let own = *b;
                            use itertools::partition;

                            let mut sec = vec!();

                            // partition to move ones inside the bound to beginning
                            if *st < ps.len() - 1 {
                                let e = partition(&mut ps[*st..], |p|{ own.inside(p) }); 
                                if e > *st {
                                    // nab the items needed for this
                                    for p in &ps[*st..e] {
                                        sec.push(*p);
                                    }
                                }
                                *st = e;
                            }
                            // make the t
                            let qtree : Self = Self::make_inner(sec, own, itr+1);

                            // move the starting point for the next partition on
                            

                            // yield the t
                            Some(qtree)
                        }).collect();
                    output.into_inner().unwrap()
                };

            Quadtree(Box::new(QuadtreeInner::Branch { split : split, children : kids }), bound)
        }
    }

    fn closest_point_inner<'a, D:Dist<V>>(&'a self, target: &V, current_best : &mut Option<(D::Output, &'a V)>) {
        let Quadtree(ref q, ref bound) = *self;

        let closest_bounded = bound.clamp(target);
       
        match *current_best {
            Some((ref d, _)) if D::dist_between(target, &closest_bounded) > *d => (),
            _ =>
                match **q {
                    QuadtreeInner::Branch{ref children, ref split} => {
                        // pick search order & enumerate
                        for idx in split.pick_search_order(target).into_iter() {
                            let child = &children[*idx];
                            child.closest_point_inner::<D>(target, current_best);
                        }
                    },
                    QuadtreeInner::Leaf{ref points, phantom : _} => {
                        for p in points {
                            let dist = D::dist_between(target, p);
                            match *current_best {
                                Some((ref d, _)) if dist > *d => (),
                                _ => *current_best = Some((dist, p))
                            }
                        }
                    },
                }
        }
    }

}

use std::iter::{Empty, Filter, FlatMap};

struct QIter<'a, V : 'a> ( Box<Iterator<Item=&'a V> + 'a> );


impl<'a, V> Iterator for QIter<'a, V> {
     type Item = &'a V;

     fn next(&mut self) -> Option<Self::Item>
     {
         self.0.next()
     }
}

impl<'a, T : Vectorable, V : Vector2<T> + 'a, B : BSPBound<V>> SpaceTree<'a, V> for Quadtree<B,V,T> {

    type I = QIter<'a, V>;

    fn new<Pts : IntoIterator<Item=V>>(pts:Pts) -> Self {
        
        let pts_clone : Vec<_> = pts.into_iter().collect();

        let bound = B::new(pts_clone.iter().cloned());

        Quadtree::make_inner(pts_clone, bound, 0)
    }

    fn closest_point<'b : 'a, D:Dist<V>>(&'a self, target : &'b V) -> Option<(D::Output, &'a V)> {
        let mut opt = Option::None;
        self.closest_point_inner::<D>(target, &mut opt);
        opt
    }

    fn points_within_dist<'b : 'a, D:Dist<V>>(&'a self, target:&'b V, dist:&'b D::Output) -> Self::I {
        let Quadtree(ref q, ref bound) = *self;

        let closest_bounded = bound.clamp(target);

        let b : Box<Iterator<Item=&'a V> + 'a> = 
        {
            if D::dist_between(&closest_bounded, target) > *dist {
                use std::iter::empty;
                Box::new(empty::<&'a V>())
            }
            else {
                match **q {
                    QuadtreeInner::Branch{ref children, ..} => {
                        Box::new(children.iter().flat_map(move |q|{ q.points_within_dist::<D>(target, dist)}))
                    },
                    QuadtreeInner::Leaf{ref points, ..} => {
                        Box::new(points.iter().filter(move |p|{ D::dist_between(p, target) < *dist}))
                    }
                }
            }

        };
        QIter(b)
    }
}

mod tests {

    extern crate rand;
    use rand::Rng;
    use rand::distributions::{IndependentSample, Range};

    use quadtree;

    type SimpleVec2<T> = [T;2];

    use vec::Vector2;
    use vec::Vectorable;

    use std::f32;    

    impl Vectorable for f32 {
        fn sqrt(&self) -> f32 { f32::sqrt(*self) }
        fn max_value() -> f32 { f32::MAX }
        fn min_value() -> f32 { -f32::MAX }
    }

    impl<T> Vector2<T> for SimpleVec2<T> where T : Vectorable
    {
        fn x(&self) -> T {
            self[0]
        }

        fn y(&self) -> T {
            self[1]
        }

        fn new(x:T, y:T) -> Self {
            [x,y]
        }
    }

    use std::marker::PhantomData;
    
    #[derive(Debug)]
    pub struct Halfway<T : Vectorable, V : Vector2<T>> (PhantomData<T>,V);

    impl<T : Vectorable, V : Vector2<T>> quadtree::BSPBoundSplit<V> for Halfway<T, V> {
        fn pick_search_order(&self, pt: &V) -> [usize; 4] {
            //work out the closest sub-quad and check that first - will hopefully reduce 'best_sq_dist' so we don't check further
            let h = self.1;
            let x_half = pt.x() > h.x(); // true if in the second half of x
            let y_half = pt.y() > h.y(); // true if in the second half of y
            let off_idx = (x_half as usize) + ((y_half as usize) << 1); // graphical representation below...
            /*
            ---------------
            |      |      |
            |  00  |  01  |
            |      |      |
            ---------------
            |      |      |
            |  10  |  11  |
            |      |      |
            ---------------
            */
            let modulate = |i| { (i + off_idx) % 4 };
            [
                modulate(0),
                modulate(1),
                modulate(2),
                modulate(3)
            ]
        }
    }

    type SimpleBound<T, V> = (PhantomData<T>, [V;2]);

    impl<T : Vectorable, V : Vector2<T>> quadtree::BSPBound<V> for SimpleBound<T, V> {

        type S = Halfway<T, V>;

        fn new<P : Iterator<Item=V>>(pts:P) -> Self {
            let mut maxx = T::min_value();
            let mut maxy = T::min_value();
            let mut minx = T::max_value();
            let mut miny = T::max_value();

            fn max<T : Vectorable>(x:T, y:T) -> T { if x > y {x} else {y} };
            fn min<T : Vectorable>(x:T, y:T) -> T { if x < y {x} else {y} };

            for p in pts.into_iter() {
                maxx = max::<T>(maxx, p.x());
                maxy = max::<T>(maxy, p.y());
                minx = min::<T>(minx, p.x());
                miny = min::<T>(miny, p.y());
            }

            (PhantomData, [V::new(minx, miny), V::new(maxx, maxy)])
        }

        fn clamp(&self, p : &V) -> V {
            let min = (self.1)[0];
            let max = (self.1)[1];

            let y : T = 
                {
                    if p.y() < min.y() {
                        min.y()
                    }
                    else if p.y() > max.y() {
                        max.y()
                    }
                    else{
                        p.y()
                    }
                };

            let x = 
                {
                    if p.x() < min.x() {
                        min.x()
                    }
                    else if p.x() > max.x() {
                        max.x()
                    }
                    else{
                        p.x()
                    }
                };

            V::new(x,y)
        }

        fn inside(&self, pt : &V) -> bool {
            let z = T::zero();

            let min = (self.1)[0];
            let max = (self.1)[1];

            (pt.x() - min.x()) >= z && (pt.x() - max.x()) <= z && (pt.y() - min.y()) >= (-z) && (pt.y() - max.y()) <= z
        }

        fn carve(&self) -> (Self::S,[Self;4])
        {
            let min = &self.1[0];
            let max = &self.1[1];

            let hx = (min.x() + max.x()) / T::two();
            let hy = (min.y() + max.y()) / T::two();

            let split = Halfway(PhantomData, V::new(hx, hy));

            let bounds = 
                [
                    (PhantomData, [*min, V::new(hx, hy)]),
                    (PhantomData, [V::new(hx, min.y()), V::new(max.x(), hy)]),
                    (PhantomData, [V::new(min.x(), hy), V::new(hx, max.y())]),
                    (PhantomData, [V::new(hx, hy), *max]),
                ];

            (split, bounds)
        }
    }

    fn random_pts(size : usize) -> Vec<SimpleVec2<f32>> {
        let mut rng = rand::thread_rng();
        let between = Range::new(-1., 1.);

        let mut v = vec![];

        for _ in 0..size {
            let x = between.ind_sample(&mut rng);
            let y = between.ind_sample(&mut rng);
            
            v.push([x,y]);
        };

        v
    }

    fn closest_point_brute_force<'a>(pts:&'a Vec<SimpleVec2<f32>>, p:&SimpleVec2<f32>) -> (f32, &'a SimpleVec2<f32>)
    {
        let mut best = &pts[0];
        let mut best_dist_sq = best.dist_sq(p);

        for pt in pts {
            let d = pt.dist_sq(p);
            if d < best_dist_sq {
                best_dist_sq = d;
                best = pt;
            }
        }

        (best_dist_sq.sqrt(), best)
    }
    
    #[derive(Debug)]
    struct SimpleVecDist;

    use spacetree::Dist;
    impl Dist<SimpleVec2<f32>> for SimpleVecDist {
        type Output = f32;
        fn dist_between(s1 : &SimpleVec2<f32>, s2 : &SimpleVec2<f32>) -> f32 {
            return s1.dist_sq(s2).sqrt();
        }
    }

    #[test]
    fn closest_point_returns_same_as_brute_force() {
        let pts = random_pts(10);

        let bf_clone = pts.clone();

        use quadtree::Quadtree;
        use spacetree::SpaceTree;

        let q = Quadtree::< SimpleBound<f32, SimpleVec2<f32> >, SimpleVec2<f32>, f32>::new(pts);

        for test_pt in random_pts(10) {
            let closest_bf = Option::Some(closest_point_brute_force(&bf_clone, &test_pt));
            let closest_q = q.closest_point::<SimpleVecDist>(&test_pt);

            assert_eq!(closest_bf, closest_q);
        }

        ()
    }

    #[test]
    fn points_within_returns_same_as_brute_force() {
        let pts = random_pts(10);

        let bf_clone = pts.clone();

        use quadtree::Quadtree;
        use spacetree::SpaceTree;

        let q = Quadtree::< SimpleBound<f32,SimpleVec2<f32> >, SimpleVec2<f32>, f32>::new(pts);

        let d = 1.0;

        for test_pt in random_pts(10) {

            use std::cmp::Ordering;

            let mut pts_within : Vec<&SimpleVec2<f32>> = bf_clone.iter().filter(|p|{ p.dist_sq(&test_pt).sqrt() < d}).collect();
            pts_within.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Less));
            let mut pts_within_test : Vec<&SimpleVec2<f32>> = q.points_within_dist::<SimpleVecDist>(&test_pt, &d).collect();
            pts_within_test.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Less));

            assert_eq!(pts_within, pts_within_test);
        }

        ()
    }

}