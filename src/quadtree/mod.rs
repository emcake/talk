
#[derive(Debug)]
#[derive(Copy)]
#[derive(Clone)]
#[derive(PartialEq)]
pub struct Vec2 {
    pub x : f32,
    pub y : f32
}

impl Vec2 {
    pub fn dist_sq(&self, pt:&Vec2) -> f32 {
        let dx = self.x - pt.x;
        let dy = self.y - pt.y;
        dx*dx + dy*dy
    }

    pub fn dist(&self, pt:&Vec2) -> f32 {
        self.dist_sq(pt).sqrt()
    }
}

#[derive(Debug)]
#[derive(Copy)]
#[derive(Clone)]
#[derive(PartialEq)]
pub struct Bound2 {
    pub min : Vec2,
    pub max : Vec2
}

impl Bound2 {

    pub fn new(x1:f32, x2:f32, y1:f32, y2:f32) -> Self {
        Bound2 { 
            min : Vec2 { x : x1, y : y1 },
            max : Vec2 { x : x2, y : y2 }
        }
    }

    pub fn split(&self, hx:f32, hy:f32) -> [Bound2;4] {
        [
            Bound2::new(self.min.x, hx, self.min.y, hy),
            Bound2::new(hx, self.max.x, self.min.y, hy),
            Bound2::new(self.min.x, hx, hy, self.max.y),
            Bound2::new(hx, self.max.x, hy, self.max.y)
        ]
    }
}

pub enum Quadtree {
    Branch {
        bound : Bound2,
        children : [Box<Quadtree>; 4]
    },
    Leaf {
        bound : Bound2,
        points : Vec<Vec2>
    }
}

mod bound {
    use quadtree::Bound2;
    use quadtree::Vec2;

    pub fn within_dist_of_bound(b : &Bound2, pt : &Vec2, rad : &f32) -> bool {
        (pt.x - b.min.x) >= -rad && (pt.x - b.max.x) <= *rad && (pt.y - b.min.y) >= -rad && (pt.y - b.max.y) <= *rad
    }

    #[test]
    fn test_within_dist_of_bound__inside() {
        let b = Bound2 {
            min : Vec2 { x : 0.0, y : 0.0},
            max : Vec2 { x : 1.0, y : 1.0}
        };

        let p1 = Vec2 { x : 0.5, y : 0.5 };

        assert_eq!(within_dist_of_bound(&b, &p1, &0.0), true);
    }

    #[test]
    fn test_within_dist_of_bound_outside() {
        let b = Bound2 {
            min : Vec2 { x : 0.0, y : 0.0},
            max : Vec2 { x : 1.0, y : 1.0}
        };

        let p1 = Vec2 { x : -0.1, y : -0.1 };

        assert_eq!(within_dist_of_bound(&b, &p1, &0.0), false);
    }

    #[test]
    fn test_within_dist_of_bound_inside_with_rad() {
        let b = Bound2 {
            min : Vec2 { x : 0.0, y : 0.0},
            max : Vec2 { x : 1.0, y : 1.0}
        };

        let p1 = Vec2 { x : -0.1, y : -0.1 };

        assert_eq!(within_dist_of_bound(&b, &p1, &0.2), true);
    }

    #[test]
    fn test_within_dist_of_bound_inside_edge() {
        let b = Bound2 {
            min : Vec2 { x : 0.0, y : 0.0},
            max : Vec2 { x : 1.0, y : 1.0}
        };

        let p1 = Vec2 { x : -0.1, y : -0.1 };

        assert_eq!(within_dist_of_bound(&b, &p1, &0.1), true);
    }

}

impl Quadtree {
    fn bound(&self) -> Bound2 {
        match *self {
            Quadtree::Branch { bound : b, children : _ } => b,
            Quadtree::Leaf {bound : b, points : _ } => b,
        }
    }

    
    pub fn closest_point_inner(&self, pt : &Vec2, mut best_sq_dist : &mut f32, mut best : &mut Vec2) -> () {
        let b = self.bound();

        let rad = best_sq_dist.sqrt();

        if bound::within_dist_of_bound(&b, &pt, &rad) {
            println!("INB: {:?} {:?}",*pt, b);
            
            match self {
                &Quadtree::Branch { bound : _, ref children } => {
                    //work out the closest sub-quad and check that first - will hopefully reduce 'best_sq_dist' so we don't check further
                    let x_half = pt.x > (b.min.x + b.min.y) / 2.0; // true if in the second half of x
                    let y_half = pt.y > (b.min.y + b.max.y) / 2.0; // true if in the second half of y
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
                        println!("CHECK_SUBQUAD_{:?}: {:?}", off_idx, children[i].bound());
                        children[i].closest_point_inner(pt, &mut best_sq_dist, &mut best);
                    }
                },
                &Quadtree::Leaf { bound : _, ref points } => {
                    for candidate in points.iter() {
                        let sq_dist = candidate.dist_sq(pt);
                        if sq_dist < *best_sq_dist {
                            *best_sq_dist = sq_dist;
                            *best = *candidate;
                            println!("NEW_BEST: {:?} @ {:?}", best_sq_dist.sqrt(), best);
                        }
                    }
                }
            }
        }
        else
        {
            println!("OOB: {:?} {:?}",*pt, b);
        }
    }

    pub fn closest_point(&self, pt : &Vec2) -> (f32,Vec2) {
        use std::f32; // you can put use declarations anywhere

        let mut best_sq_dist = f32::MAX;

        let mut best_point = Vec2 { x : f32::NAN, y : f32 ::NAN };

        self.closest_point_inner(pt, &mut best_sq_dist, &mut best_point);

        (best_sq_dist.sqrt(), best_point)
    }

    pub fn points_close_to(&self, pt: &Vec2, rad:f32, pts: &mut Vec<Vec2>) -> ()
    {
        let b = self.bound();

        if (pt.x - b.min.x) < rad && (b.max.x - pt.x) > -rad {
            match self {
                &Quadtree::Branch { bound : _, ref children } => {
                    for subtree in children.iter() {
                        subtree.points_close_to(pt, rad, pts);
                    }
                },
                &Quadtree::Leaf { bound : _, ref points } => {
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

    fn lerp (x:f32, y:f32, amt:f32) -> f32 {
        x + ((y-x)*amt)
    }
        

    fn make_inner(pts : &Vec<Vec2>, bound : Bound2, itr : i32) -> Box<Quadtree> {
        if pts.len() <= 1 || itr >= 10 {
            Box::new(Quadtree::Leaf {bound : bound, points : pts.clone()})
        }
        else {
            
            let hx = Quadtree::lerp(bound.min.x, bound.max.x, 0.5);
            let hy = Quadtree::lerp(bound.min.y, bound.max.y, 0.5);

            let mut top_left : Vec<Vec2> = vec![];
            let mut top_right : Vec<Vec2> = vec![];
            let mut bot_left : Vec<Vec2> = vec![];
            let mut bot_right : Vec<Vec2> = vec![];

            for p in pts.iter() {
                match (p.x < hx, p.y < hy) {
                    (true, true) => top_left.push(*p),
                    (false, true) => top_right.push(*p),
                    (true, false) => bot_left.push(*p),
                    (false, false) => bot_right.push(*p)
                }              
            };

            let bounds = bound.split(hx, hy);

            let kids = [
                Quadtree::make_inner(&top_left, bounds[0], itr + 1),
                Quadtree::make_inner(&top_right, bounds[1], itr + 1),
                Quadtree::make_inner(&bot_left, bounds[2], itr + 1),
                Quadtree::make_inner(&bot_right, bounds[3], itr + 1),
            ];

            Box::new(Quadtree::Branch { bound : bound, children : kids })
        }
    }

    pub fn make(pts : &Vec<Vec2>) -> Box<Quadtree> {
        let mut maxx = 0. / 0.;
        let mut maxy = 0. / 0.;
        let mut minx = 1. / 0.;
        let mut miny = 1. / 0.;

        for p in pts.iter() {
            maxx = f32::max(maxx, p.x);
            maxy = f32::max(maxy, p.y);
            minx = f32::min(minx, p.x);
            miny = f32::min(miny, p.y);
        }

        Quadtree::make_inner(pts, Bound2 { min : Vec2 { x : minx, y : miny }, max : Vec2 { x : maxx, y : maxy } }, 0)
    }

}

mod tests {

    extern crate rand;
    use rand::Rng;
    use rand::distributions::{IndependentSample, Range};

    use quadtree;
    use quadtree::Vec2;

    fn random_pts(size : usize) -> Vec<Vec2> {
        let mut rng = rand::thread_rng();
        let between = Range::new(-1., 1.);

        let mut v = vec![];

        for _ in 0..size {
            let x = between.ind_sample(&mut rng);
            let y = between.ind_sample(&mut rng);
            
            v.push(Vec2{x : x, y : y});
        };

        v
    }

    fn closest_point_brute_force(pts:&Vec<Vec2>, p:&Vec2) -> (f32, Vec2)
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

        let q = Quadtree::make(&pts);

        for test_pt in random_pts(10) {
            let closest_bf = closest_point_brute_force(&pts, &test_pt);
            let closest_q = q.closest_point(&test_pt);

            assert_eq!(closest_bf, closest_q);
        }

        ()
    }

}