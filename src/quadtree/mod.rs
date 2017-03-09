
#[derive(Debug)]
#[derive(Copy)]
#[derive(Clone)]
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
pub struct Bound2 {
    pub min : Vec2,
    pub max : Vec2
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

impl Quadtree {
    fn bound(&self) -> Bound2 {
        match *self {
            Quadtree::Branch { bound : b, children : _ } => b,
            Quadtree::Leaf {bound : b, points : _ } => b,
        }
    }

    fn closest_point_inner(&self, pt : &Vec2, mut best_sq_dist : &mut f32, mut best : &mut Vec2) -> () {
        let b = self.bound();

        let rad = best_sq_dist.sqrt();

        if (pt.x - b.min.x) < rad && (b.max.x - pt.x) > -rad {
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
                        children[i].closest_point_inner(pt, &mut best_sq_dist, &mut best);
                    }
                },
                &Quadtree::Leaf { bound : _, ref points } => {
                    for candidate in points.iter() {
                        let sq_dist = candidate.dist_sq(pt);
                        if sq_dist < *best_sq_dist {
                            *best_sq_dist = sq_dist;
                            *best = *candidate;
                        }
                    }
                }
            }
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

    fn make_inner(pts : &Vec<Vec2>, bound : Bound2, itr : i32) -> Box<Quadtree> {
        if pts.len() <= 1 || itr >= 10 {
            Box::new(Quadtree::Leaf {bound : bound, points : pts.clone()})
        }
        else {
            let split_x = (bound.min.x + bound.max.x) / 2.;
            let split_y = (bound.min.y + bound.max.y) / 2.;

            let mut top_left : Vec<Vec2> = vec![];
            let mut top_right : Vec<Vec2> = vec![];
            let mut bot_left : Vec<Vec2> = vec![];
            let mut bot_right : Vec<Vec2> = vec![];

            for p in pts.iter() {
                match (p.x < split_x, p.y < split_y) {
                    (true, true) => top_left.push(*p),
                    (false, true) => top_right.push(*p),
                    (true, false) => bot_left.push(*p),
                    (false, false) => bot_right.push(*p)
                }              
            };

            let kids = [
                Quadtree::make_inner(&top_left, Bound2 { min : bound.min, max : Vec2 { x : split_x, y : split_y } }, itr + 1),
                Quadtree::make_inner(&top_right, Bound2 { min : Vec2{ x : split_x, y : bound.min.y }, max : Vec2 { x : bound.max.x, y : split_y } }, itr + 1),
                Quadtree::make_inner(&bot_left, Bound2 { min : Vec2{ x : bound.min.x, y : split_y }, max : Vec2 { x : split_x, y : bound.max.y } }, itr + 1),
                Quadtree::make_inner(&top_left, Bound2 { min : Vec2{ x : split_x, y : split_y }, max : bound.max }, itr + 1),
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

        Quadtree::make_inner(pts, Bound2 { min : Vec2 { x : minx, y : miny }, max : Vec2 { x : minx, y : miny } }, 0)
    }

}

#[cfg(test)]
mod tests {

    extern crate rand;
    use rand::Rng;
    use rand::distributions::{IndependentSample, Range};

    use quadtree::test::Bencher;

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

}