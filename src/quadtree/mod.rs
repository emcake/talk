extern crate arrayvec;

use vec::Vector2;
use vec::Bound2;
use vec::Vectorable;

use std::marker::PhantomData;
    
#[derive(Debug)]
pub enum Quadtree<B,V,T> where B : Bound2<V,T>, V : Vector2<T>, T : Vectorable {
    Branch {
        bound : B,
        children : [Box<Quadtree<B,V,T>>; 4]
    },
    Leaf {
        bound : B,
        points : Vec<V>,
        phantom : PhantomData<T>
    }
}

impl<B,V,T> Quadtree<B,V,T> where B : Bound2<V,T>, V : Vector2<T>, T : Vectorable {
    fn bound(&self) -> &B {
        match self {
            &Quadtree::Branch { bound : ref b, children : _ } => b,
            &Quadtree::Leaf {bound : ref b, points : _ , phantom: _} => b,
        }
    }

    
    pub fn closest_point_inner(&self, pt : &V, mut best_sq_dist : &mut T, mut best : &mut Option<V>) -> () {
        let b = self.bound();

        let rad = best_sq_dist.sqrt();

        if b.point_within(pt, rad) {
            //println!("INB: {:?} {:?}",*pt, b);
            
            match self {
                &Quadtree::Branch { bound : _, ref children } => {
                    //work out the closest sub-quad and check that first - will hopefully reduce 'best_sq_dist' so we don't check further
                    let x_half = pt.x() > (b.min().x() + b.min().y()) / T::two(); // true if in the second half of x
                    let y_half = pt.y() > (b.min().y() + b.max().y()) / T::two(); // true if in the second half of y
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
                    for i in 0..4 {
                        let i = (i + off_idx) % 4;
                        //println!("CHECK_SUBQUAD_{:?}: {:?}", off_idx, children[i].bound());
                        children[i].closest_point_inner(pt, &mut best_sq_dist, &mut best);
                    }
                },
                &Quadtree::Leaf { bound : _, ref points, phantom : _ } => {
                    for candidate in points.iter() {
                        let sq_dist = candidate.dist_sq(pt);
                        if sq_dist < *best_sq_dist {
                            *best_sq_dist = sq_dist;
                            *best = Option::Some(candidate.clone());
                            //println!("NEW_BEST: {:?} @ {:?}", best_sq_dist.sqrt(), best);
                        }
                    }
                }
            }
        }
        else
        {
            //println!("OOB: {:?} {:?}",*pt, b);
        }
    }

    pub fn closest_point(&self, pt : &V) -> Option<(T,V)> {
        use std::f32; // you can put use declarations anywhere

        let mut best_sq_dist = T::max_value();

        let mut best_point = None;

        self.closest_point_inner(pt, &mut best_sq_dist, &mut best_point);

        best_point.map(|p| -> (T,V) { (best_sq_dist.sqrt(), p) })
    }

    pub fn points_close_to(&self, pt: &V, rad:T, pts: &mut Vec<V>) -> ()
    {
        let b = self.bound();

        if (pt.x() - b.min().x()) < rad && (b.max().x() - pt.x()) > -rad {
            match self {
                &Quadtree::Branch { bound : _, ref children } => {
                    for subtree in children.iter() {
                        subtree.points_close_to(pt, rad, pts);
                    }
                },
                &Quadtree::Leaf { bound : _, ref points, phantom : _ } => {
                    let rad_sq = rad*rad;
                    for candidate in points.iter() {
                        if candidate.dist_sq(pt) <= rad_sq {
                            pts.push(*candidate);
                        }
                    }
                }
            }
        }
    }

    
        

    fn make_inner(pts : &Vec<V>, bound : B, itr : i32) -> Box<Quadtree<B,V,T>> {
        if pts.len() <= 1 || itr >= 10 {
            Box::new(Quadtree::Leaf {bound : bound, points : pts.clone(), phantom : PhantomData})
        }
        else {

            fn half<T : Vectorable> (x:T, y:T) -> T {
                x + ((y-x)/ T::two())
            }
            
            let hx = half(bound.min().x(), bound.max().x());
            let hy = half(bound.min().y(), bound.max().y());

            let mut vecs = [vec![], vec![], vec![], vec![]];

            for p in pts.iter() {
                match (p.x() < hx, p.y() < hy) {
                    (true, true) => vecs[0].push(*p),
                    (false, true) => vecs[1].push(*p),
                    (true, false) => vecs[2].push(*p),
                    (false, false) => vecs[3].push(*p)
                }              
            };

            let bounds = bound.split(hx, hy);

            use self::arrayvec::ArrayVec;

            let kids : [Box<Quadtree<B,V,T>>;4] = 
                {
                    let output : ArrayVec<[Box<Quadtree<B,V,T>>;4]> = 
                        bounds.into_iter().zip(vecs.into_iter()).map(|(&b,v)| -> Box<Quadtree<B,V,T>> {
                            Quadtree::make_inner(v, b, itr+1)
                        }).collect();
                    output.into_inner().unwrap()
                };

            Box::new(Quadtree::Branch { bound : bound, children : kids })
        }
    }

    pub fn make(pts : &Vec<V>) -> Box<Quadtree<B,V,T>> {
        let mut maxx = T::min_value();
        let mut maxy = T::min_value();
        let mut minx = T::max_value();
        let mut miny = T::max_value();

        fn max<T : Vectorable>(x:T, y:T) -> T { if x > y {x} else {y} };
        fn min<T : Vectorable>(x:T, y:T) -> T { if x < y {x} else {y} };

        for p in pts.iter() {
            maxx = max::<T>(maxx, p.x());
            maxy = max::<T>(maxy, p.y());
            minx = min::<T>(minx, p.x());
            miny = min::<T>(miny, p.y());
        }

        Quadtree::make_inner(pts, B::new(V::new(minx, miny), V::new(maxx, maxy)), 0)
    }

    
}

mod tests {

    extern crate rand;
    use rand::Rng;
    use rand::distributions::{IndependentSample, Range};

    use quadtree;

    type SimpleVec2<T> = [T;2];

    use vec::Vector2;
    use vec::Bound2;
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

    type SimpleBound<V> = [V;2];

    impl<V : Vector2<T>,T : Vectorable> Bound2<V,T> for SimpleBound<V> {
        fn min(&self) -> &V { &self[0] }
        fn max(&self) -> &V { &self[1] }

        fn new(mn:V, mx:V) -> Self {
            [mn,mx]
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

    fn closest_point_brute_force(pts:&Vec<SimpleVec2<f32>>, p:&SimpleVec2<f32>) -> (f32, SimpleVec2<f32>)
    {
        let mut best = pts[0];
        let mut best_dist_sq = best.dist_sq(p);

        for pt in pts {
            let d = pt.dist_sq(p);
            if d < best_dist_sq {
                best_dist_sq = d;
                best = *pt;
            }
        }

        (best_dist_sq.sqrt(), best)
    }

    #[test]
    fn closest_point_returns_same_as_brute_force() {
        let pts = random_pts(100);

        use quadtree::Quadtree;

        let q = Quadtree::< SimpleBound<SimpleVec2<f32> >, SimpleVec2<f32>, f32>::make(&pts);

        for test_pt in random_pts(10) {
            let closest_bf = Option::Some(closest_point_brute_force(&pts, &test_pt));
            let closest_q = q.closest_point(&test_pt);
            assert!(false);

            assert_eq!(closest_bf, closest_q);
        }

        ()
    }

}