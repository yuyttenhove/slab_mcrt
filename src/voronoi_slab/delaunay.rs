use super::geometry::{orient_2d, in_circle_2d, circumcenter_2d, circumradius_2d};
use super::NeighbourDirection;
use crate::vector::Vec2;
use std::collections::VecDeque;
use std::fs;
use permutation::Permutation;
use crate::voronoi_slab::voronoi_generator::VoronoiGenerator;

pub fn random_choose<T>(option1: T, option2: T) -> T {
    if rand::random() {option1} else {option2}
}

#[derive(Debug)]
pub(super) struct DelaunayVertex {
    pub(super) x: f64,
    pub(super) y: f64,
    pub(super) x_scaled: f64,
    pub(super) y_scaled: f64,
    pub(super) triangle: i32,
    pub(super) index_in_triangle: i8,
    search_radius: f64
}

impl Default for DelaunayVertex {
    fn default() -> DelaunayVertex {
        DelaunayVertex {
            x: f64::NAN,
            y: f64::NAN,
            x_scaled: f64::NAN,
            y_scaled: f64::NAN,
            triangle: -1,
            index_in_triangle: -1,
            search_radius: f64::INFINITY
        }
    }
}

impl DelaunayVertex {
    fn update_triangle(&mut self, triangle: i32, triangle_index: i8) {
        self.triangle = triangle;
        self.index_in_triangle = triangle_index;
    }
}


#[derive(Debug)]
pub(super) struct DelaunayTriangle {
    pub(super) vertices: [i32; 3],
    pub(super) neighbours: [i32; 3],
    pub(super) index_in_neighbours: [i8; 3]
}

impl Default for DelaunayTriangle {
    fn default() -> DelaunayTriangle {
        DelaunayTriangle {
            vertices: [-1, -1, -1],
            neighbours: [-1, -1, -1],
            index_in_neighbours: [-1, -1, -1]
        }
    }
}

impl DelaunayTriangle {
    fn update_neighbours(&mut self, n0: i32, n1: i32, n2: i32, idx_in_n0: i8, idx_in_n1: i8, idx_in_n2: i8) {
        self.neighbours = [n0, n1, n2];
        self.index_in_neighbours = [idx_in_n0, idx_in_n1, idx_in_n2];
    }

    fn update_neighbour(&mut self, n: i32, idx_in_n: i8, i: i8) {
        self.neighbours[i as usize] = n;
        self.index_in_neighbours[i as usize] = idx_in_n;
    }

    pub(super) fn circumcenter(&self, triangulation: &DelaunayTriangulation) -> Vec2<f64> {
        circumcenter_2d(
            triangulation.vertices[self.vertices[0] as usize].x,
            triangulation.vertices[self.vertices[0] as usize].y,
            triangulation.vertices[self.vertices[1] as usize].x,
            triangulation.vertices[self.vertices[1] as usize].y,
            triangulation.vertices[self.vertices[2] as usize].x,
            triangulation.vertices[self.vertices[2] as usize].y
        )
    }

    fn circumradius(&self, triangulation: &DelaunayTriangulation) -> f64 {
        circumradius_2d(
            triangulation.vertices[self.vertices[0] as usize].x,
            triangulation.vertices[self.vertices[0] as usize].y,
            triangulation.vertices[self.vertices[1] as usize].x,
            triangulation.vertices[self.vertices[1] as usize].y,
            triangulation.vertices[self.vertices[2] as usize].x,
            triangulation.vertices[self.vertices[2] as usize].y,
        )
    }
}


#[derive(Debug, Default)]
pub struct DelaunayTriangulation {
    pub(super) vertices: Vec<DelaunayVertex>,
    pub(super) triangles: Vec<DelaunayTriangle>,
    pub(super) is_periodic: bool,
    pub(super) n_vertices: usize,
    pub(super) anchor: Vec2<f64>,
    pub(super) sides: Vec2<f64>,
    pub(super) neighbour_info: Vec<(usize, NeighbourDirection)>,
    side: f64,
    inverse_side: f64,
    current_triangle_idx: i32,
    current_vertex_idx: i32,
    triangles_to_check: VecDeque<i32>
}

impl DelaunayTriangulation {
    pub fn new(anchor: Vec2<f64>, sides: Vec2<f64>, vertex_size: usize, triangle_size: usize) -> DelaunayTriangulation {
        let mut triangulation = DelaunayTriangulation {
            vertices: Vec::with_capacity(vertex_size + 3),
            triangles: Vec::with_capacity(triangle_size * 2),
            ..DelaunayTriangulation::default()
        };

        /* Setup the domain and side of the triangulation box large enough so that any
        ghost particles certainly fall into the domain.
         */
        triangulation.anchor = anchor - sides;
        triangulation.sides = sides;

        triangulation.side = 6. * f64::max(sides.x(), sides.y());
        triangulation.inverse_side = 1. / triangulation.side;

        // Create vertices for the first triangle that encapsulates the entire domain
        let v0 = triangulation.new_vertex(triangulation.anchor.x(), triangulation.anchor.y());
        let v1 = triangulation.new_vertex(triangulation.anchor.x() + triangulation.side,
                                          triangulation.anchor.y());
        let v2 = triangulation.new_vertex(triangulation.anchor.x(),
                                          triangulation.anchor.y() + triangulation.side);

        // Create first large triangle and 3 dummy triangles
        let dummy0 = triangulation.new_triangle(v1, -1, v2);
        let dummy1 = triangulation.new_triangle(v2, -1, v0);
        let dummy2 = triangulation.new_triangle(v0, -1, v1);
        let first_triangle = triangulation.new_triangle(v0, v1, v2);
        // set neighbour relations: first triangles has 3 dummies as neighbours and each dummy has
        // the first triangle and the other two dummies as neighbours.
        triangulation.triangles[first_triangle as usize].update_neighbours(dummy0, dummy1, dummy2,
                                                                           1, 1, 1);
        triangulation.triangles[dummy0 as usize].update_neighbours(dummy1, first_triangle, dummy2,
                                                                   2, 0, 0);
        triangulation.triangles[dummy1 as usize].update_neighbours(dummy2, first_triangle, dummy0,
                                                                   2, 1, 0);
        triangulation.triangles[dummy2 as usize].update_neighbours(dummy0, first_triangle, dummy1,
                                                                   2, 2, 0);

        triangulation.consistency_check();

        triangulation.current_triangle_idx = first_triangle;
        triangulation
    }

    pub fn from_points(points_x: &Vec<f64>,
                       points_y: &Vec<f64>,
                       anchor: Vec2<f64>,
                       sides: Vec2<f64>,
                       make_periodic: bool) -> DelaunayTriangulation {
        assert_eq!(points_x.len(), points_y.len(), "points_x and points_y must have the same length!");
        let mut d = DelaunayTriangulation::new(anchor, sides, points_x.len(), points_y.len() * 2);
        for (i, &x) in points_x.iter().enumerate() {
            d.insert_point(x, points_y[i]);
        }
        d.n_vertices = d.vertices.len() - 3;
        if make_periodic {
            d.make_periodic();
        }
        d
    }

    pub fn from_generators(generators: &Vec<VoronoiGenerator>, anchor: Vec2<f64>, sides: Vec2<f64>, make_periodic: bool) -> DelaunayTriangulation {
        let mut d = DelaunayTriangulation::new(anchor, sides, generators.len(), 6*generators.len());
        for  gen in generators {
            d.insert_point(gen.x(), gen.y())
        }
        d.n_vertices = d.vertices.len() - 3;
        if make_periodic {
            d.make_periodic();
        }
        d
    }

    fn insert_point(&mut self, x: f64, y: f64) {
        // add vertex
        self.new_vertex(x, y);

        // Find triangle in which (x, y) is positioned
        self.current_triangle_idx = self.find_triangle_containing_current_vertex();
        let triangle = &self.triangles[self.current_triangle_idx as usize];

        // Create 3 new triangles
        let (v0, v1, v2) = (triangle.vertices[0], triangle.vertices[1], triangle.vertices[2]);
        let (n0, n1, n2) = (triangle.neighbours[0], triangle.neighbours[1], triangle.neighbours[2]);
        let (idx_in_n0, idx_in_n1, idx_in_n2) = (
            triangle.index_in_neighbours[0],
            triangle.index_in_neighbours[1],
            triangle.index_in_neighbours[2]
        );

        let triangle0 = self.new_triangle_at(v1, v2, self.current_vertex_idx, self.current_triangle_idx);
        let triangle1 = self.new_triangle(v2, v0, self.current_vertex_idx);
        let triangle2 = self.new_triangle(v0, v1, self.current_vertex_idx);

        // Update neighbours
        self.triangles[triangle0 as usize].update_neighbours(triangle1, triangle2, n0,
                                                             1, 0, idx_in_n0);
        self.triangles[triangle1 as usize].update_neighbours(triangle2, triangle0, n1,
                                                             1, 0, idx_in_n1);
        self.triangles[triangle2 as usize].update_neighbours(triangle0, triangle1, n2,
                                                             1, 0, idx_in_n2);
        self.triangles[n0 as usize].update_neighbour(triangle0, 2, idx_in_n0);
        self.triangles[n1 as usize].update_neighbour(triangle1, 2, idx_in_n1);
        self.triangles[n2 as usize].update_neighbour(triangle2, 2, idx_in_n2);

        // Add new triangles to queue to check for delaunayness and check the criterion
        self.triangles_to_check.push_back(triangle0);
        self.triangles_to_check.push_back(triangle1);
        self.triangles_to_check.push_back(triangle2);

        self.fix_delaunayness();

        self.consistency_check();
    }

    fn make_periodic(&mut self) {
        use ordered_float::OrderedFloat;
        assert!(!self.is_periodic, "Delaunay triangulation is already periodic!");
        let arg_sort_x = permutation::sort_by_key(&self.vertices[3..],
                                                  |v| OrderedFloat(v.x_scaled));
        let arg_sort_y = permutation::sort_by_key(&self.vertices[3..],
                                                  |v| OrderedFloat(v.y_scaled));
        let arg_sort_xpy = permutation::sort_by_key(&self.vertices[3..],
                                                    |v| OrderedFloat(v.x_scaled + v.y_scaled));
        let arg_sort_xmy = permutation::sort_by_key(&self.vertices[3..],
                                                    |v| OrderedFloat(v.x_scaled - v.y_scaled));

        // initial value of search radius: the average inter-particle distance for uniform distribution of particles
        let mut search_radius = self.side / f64::sqrt(self.n_vertices as f64);

        let mut old_search_radius = 0.;
        let mut n_vertices_larger_search_radius = self.n_vertices;

        while n_vertices_larger_search_radius > 0 {
            self.add_ghost_vertices(search_radius, old_search_radius, &arg_sort_x, &arg_sort_y, &arg_sort_xpy, &arg_sort_xmy);
            self.update_vertex_search_radii(search_radius, n_vertices_larger_search_radius);
            n_vertices_larger_search_radius = self.vertices[3..self.n_vertices+3].iter()
                .filter(|v| (*v).search_radius > search_radius).count();
            let new_search_radius = 1.5 * search_radius;
            old_search_radius = search_radius;
            search_radius = new_search_radius;
        }

        self.is_periodic = true;
    }

    fn add_ghost_vertices(&mut self, search_radius: f64, old_search_radius: f64,
                          arg_sort_x: &Permutation, arg_sort_y: &Permutation,
                          arg_sort_xpy: &Permutation, arg_sort_xmy: &Permutation) {
        let sides = [self.sides.x(), self.sides.y()];
        // add ghost particles in positive x direction
        self.add_ghost_vertices_along_axis(
            &arg_sort_x, |v| v.x, 1,
            old_search_radius, search_radius, |v| (v.x + sides[0], v.y), NeighbourDirection::Right
        );
        // add ghost particles in negative x direction
        self.add_ghost_vertices_along_axis(
            &arg_sort_x, |v| sides[0] - v.x, -1,
            old_search_radius, search_radius, |v| (v.x - sides[0], v.y), NeighbourDirection::Left
        );
        // add ghost particles in positive y direction
        self.add_ghost_vertices_along_axis(
            &arg_sort_y, |v| sides[1] - v.y, -1,
            old_search_radius, search_radius, |v| (v.x, 2. * sides[1] - v.y), NeighbourDirection::Up
        );
        // add ghost particles in negative y direction
        self.add_ghost_vertices_along_axis(
            &arg_sort_y, |v| v.y, 1,
            old_search_radius, search_radius, |v| (v.x, -v.y), NeighbourDirection::Down
        );
        // In order to search for all particles up to r away from the corner in the diagonal
        // directions, we must compensate with a factor 1/sqrt(2); without it we could miss some
        // particles between sqrt(2)/2*r and r away from the corner...
        let sqrt2 = f64::sqrt(2.);
        // add ghost particles in positive xpy direction
        self.add_ghost_vertices_along_axis(
            &arg_sort_xmy, |v| (v.x + sides[1]- v.y) / sqrt2, 1,
            old_search_radius, search_radius, |v| (v.x + sides[0], 2. * sides[1] - v.y), NeighbourDirection::Up
        );
        // add ghost particles in negative xpy direction
        self.add_ghost_vertices_along_axis(
            &arg_sort_xmy, |v| (sides[0] - v.x + v.y) / sqrt2, -1,
            old_search_radius, search_radius, |v| (v.x - sides[0], -v.y), NeighbourDirection::Down
        );
        // add ghost particles in positive xmy direction
        self.add_ghost_vertices_along_axis(
            &arg_sort_xpy, |v| (v.x + v.y) / sqrt2, 1,
            old_search_radius, search_radius, |v| (v.x + sides[0], -v.y), NeighbourDirection::Down
        );
        // add ghost particles in negative xmy direction
        self.add_ghost_vertices_along_axis(
            &arg_sort_xpy, |v| (sides[0] - v.x + sides[1] - v.y) / sqrt2, -1,
            old_search_radius, search_radius, |v| (v.x - sides[0], 2. * sides[1] - v.y), NeighbourDirection::Up
        );
    }

    fn add_ghost_vertices_along_axis(&mut self,
                                     ordering: &Permutation,
                                     comp_f: impl Fn(&DelaunayVertex)->f64,
                                     direction: i8,
                                     old_search_radius: f64,
                                     search_radius: f64,
                                     insert_f: impl Fn(&DelaunayVertex)->(f64, f64),
                                     neighbour_direction: NeighbourDirection) {
        let mut i: usize;
        if direction == 1 { i = 0; }
        else if direction == -1 { i = self.n_vertices - 1; }
        else { panic!("The only alowed values for direction are 1 or -1!"); }

        let mut mirrored_vertex_idx = ordering.apply_inv_idx(i) + 3;
        let mut vertex = &self.vertices[mirrored_vertex_idx];
        let mut comp_value = comp_f(vertex);
        while comp_value < search_radius {
            if old_search_radius <= comp_value {
                let (x, y) = insert_f(vertex);
                self.insert_point(x, y);
                self.neighbour_info.push((mirrored_vertex_idx, neighbour_direction))
            }
            if direction == 1 {
                if i == self.n_vertices - 1 {break;}
                i += 1;
            } else {
                if i == 0 {break;}
                i -= 1;
            }
            mirrored_vertex_idx = ordering.apply_inv_idx(i) + 3;
            vertex = &self.vertices[mirrored_vertex_idx];
            comp_value = comp_f(vertex);
        }
    }

    fn update_vertex_search_radii(&mut self, current_search_radius: f64, previous_n_vertices_larger_radius: usize) {
        let mut radii = Vec::<(usize, f64)>::with_capacity(previous_n_vertices_larger_radius);
        for (i, vertex) in self.vertices[3..self.n_vertices+3].iter().enumerate() {
            if vertex.search_radius < current_search_radius { continue; }
            let mut max_radius: f64 = 0.;
            for triangle_idx in self.get_triangle_idx_around_vertex(i + 3) {
                // any point within a circle of radius 2*current circumradius around current vertex
                // CAN violate the delaunay criterion for the current triangle. In most cases
                // however this rather is unlikely...
                max_radius = max_radius.max(2. * self.triangles[triangle_idx].circumradius(self));
            }
            radii.push((i+3, max_radius));
        }
        for (i, radius) in radii {
            self.vertices[i].search_radius = radius;
        }
    }

    pub fn get_triangle_idx_around_vertex(&self, vertex_idx: usize) -> Vec::<usize> {
        let vertex = &self.vertices[vertex_idx];
        let start_triangle_idx_in_d = vertex.triangle;
        let mut triangle_indices = vec![start_triangle_idx_in_d as usize];
        let mut current_triangle_idx_in_d = start_triangle_idx_in_d;
        let mut idx_in_current_triangle = vertex.index_in_triangle;
        let mut current_triangle = &self.triangles[current_triangle_idx_in_d as usize];
        let mut next_triangle_idx_in_current_triangle = ((idx_in_current_triangle + 1) % 3) as usize;
        let mut next_triangle_idx_in_d = current_triangle.neighbours[next_triangle_idx_in_current_triangle];
        let mut current_triangle_idx_in_next_triangle = current_triangle.index_in_neighbours[next_triangle_idx_in_current_triangle];
        let mut idx_in_next_triangle = (current_triangle_idx_in_next_triangle + 1) % 3;

        while next_triangle_idx_in_d != start_triangle_idx_in_d {
            assert_eq!(
                current_triangle.vertices[idx_in_current_triangle as usize],
                self.triangles[next_triangle_idx_in_d as usize].vertices[idx_in_next_triangle as usize],
                "Vertex in next triangle is not the same as the current vertex!"
            );
            triangle_indices.push(next_triangle_idx_in_d as usize);

            current_triangle_idx_in_d = next_triangle_idx_in_d;
            idx_in_current_triangle = idx_in_next_triangle;
            current_triangle = &self.triangles[current_triangle_idx_in_d as usize];
            next_triangle_idx_in_current_triangle = ((idx_in_current_triangle + 1) % 3) as usize;
            next_triangle_idx_in_d = current_triangle.neighbours[next_triangle_idx_in_current_triangle];
            current_triangle_idx_in_next_triangle = current_triangle.index_in_neighbours[next_triangle_idx_in_current_triangle];
            idx_in_next_triangle = (current_triangle_idx_in_next_triangle + 1) % 3;
        }
        triangle_indices
    }

    pub fn is_connected_to_non_dummy_non_ghost_vertex(&self, vertex_idx: usize) -> bool {
        if vertex_idx < self.n_vertices + 3 && vertex_idx > 2 {
            return true;
        }
        let vertex = &self.vertices[vertex_idx];
        let mut idx_in_current_triangle = vertex.index_in_triangle;
        for current_triangle_idx_in_d in self.get_triangle_idx_around_vertex(vertex_idx) {
            let current_triangle = &self.triangles[current_triangle_idx_in_d];
            let next_triangle_idx_in_current_triangle = ((idx_in_current_triangle + 1) % 3) as usize;
            let neighbouring_vertex_idx = current_triangle.vertices[((idx_in_current_triangle + 2) % 3) as usize] as usize;
            if neighbouring_vertex_idx < self.n_vertices + 3 && neighbouring_vertex_idx > 2 {
                return true;
            }
            let current_triangle_idx_in_next_triangle = current_triangle.index_in_neighbours[next_triangle_idx_in_current_triangle];
            idx_in_current_triangle = (current_triangle_idx_in_next_triangle + 1) % 3;
        }
        false
    }

    fn new_vertex(&mut self, x: f64, y: f64) -> i32 {
        // TODO possibly manage the size of self.vertices more intelligently.
        let (x_scaled, y_scaled) = (x * self.inverse_side, y* self.inverse_side);
        self.current_vertex_idx = self.vertices.len() as i32;
        self.vertices.push(DelaunayVertex { x, y, x_scaled, y_scaled, ..DelaunayVertex::default()});

        self.current_vertex_idx
    }

    fn new_triangle(&mut self, v0: i32, v1: i32, v2: i32) -> i32 {
        self.new_triangle_at(v0, v1, v2, -1)
    }

    fn new_triangle_at(&mut self, v0: i32, v1: i32, v2: i32, mut at: i32) -> i32 {
        // TODO possibly manage the size of self.triangles more intelligently.
        if at < 0 {
            at = self.triangles.len() as i32;
            self.triangles.push(DelaunayTriangle {
                vertices: [v0, v1, v2],
                ..DelaunayTriangle::default()
            });
        } else {
            self.triangles[at as usize] = DelaunayTriangle {
                vertices: [v0, v1, v2],
                ..DelaunayTriangle::default()
            };
        }

        if v0 >= 0 {
            self.vertices[v0 as usize].update_triangle(at, 0);
        }
        if v1 >= 0 {
            self.vertices[v1 as usize].update_triangle(at, 1);
        }
        if v2 >= 0 {
            self.vertices[v2 as usize].update_triangle(at, 2);
        }

        at
    }

    fn find_triangle_containing_current_vertex(&self) -> i32 {
        let vertex = &self.vertices[self.current_vertex_idx as usize];
        let mut v0: &DelaunayVertex;
        let mut v1: &DelaunayVertex;
        let mut v2: &DelaunayVertex;
        let mut current_triangle: &DelaunayTriangle;
        let mut current_triangle_idx = self.current_triangle_idx;
        let mut found = false;
        let mut test0: f64;
        let mut test1: f64;
        let mut test2: f64;

        while !found {
            current_triangle = &self.triangles[current_triangle_idx as usize];
            v0 = &self.vertices[current_triangle.vertices[0] as usize];
            v1 = &self.vertices[current_triangle.vertices[1] as usize];
            v2 = &self.vertices[current_triangle.vertices[2] as usize];

            test2 = orient_2d(v0.x_scaled, v0.y_scaled, v1.x_scaled, v1.y_scaled, vertex.x_scaled, vertex.y_scaled);
            test0 = orient_2d(v1.x_scaled, v1.y_scaled, v2.x_scaled, v2.y_scaled, vertex.x_scaled, vertex.y_scaled);
            test1 = orient_2d(v2.x_scaled, v2.y_scaled, v0.x_scaled, v0.y_scaled, vertex.x_scaled, vertex.y_scaled);

            if (test0 > 0.) && (test1 > 0.) && (test2 > 0.) {
                found = true;
            } else if (test0 <= 0.) && (test1 <= 0.) {
                current_triangle_idx = random_choose(current_triangle.neighbours[0],
                                                     current_triangle.neighbours[1]);
            } else if (test1 <= 0.) && (test2 <= 0.) {
                current_triangle_idx = random_choose(current_triangle.neighbours[1],
                                                     current_triangle.neighbours[2]);
            } else if (test0 <= 0.) && (test2 <= 0.) {
                current_triangle_idx = random_choose(current_triangle.neighbours[0],
                                                     current_triangle.neighbours[2]);
            } else if test0 <= 0. {
                current_triangle_idx = current_triangle.neighbours[0];
            } else if test1 <= 0.{
                current_triangle_idx = current_triangle.neighbours[1];
            } else if test2 <= 0.{
                current_triangle_idx = current_triangle.neighbours[2];
            } else {
                panic!("Error! This scenario should not be possible?")
            }

            assert!(current_triangle_idx > 2, "Ended up with dummy triangle!");
        }
        current_triangle_idx
    }

    fn fix_delaunayness(&mut self) {
        while !self.triangles_to_check.is_empty() {
            let triangle_to_fix_idx = self.triangles_to_check.pop_front().unwrap();
            self.fix_triangle(triangle_to_fix_idx);
        }
    }

    fn fix_triangle(&mut self, triangle_idx: i32) {
        let triangle = &self.triangles[triangle_idx as usize];
        let neighbour_idx = triangle.neighbours[2];
        if neighbour_idx < 3 {
            return
        }

        let a = &self.vertices[triangle.vertices[0] as usize];
        let b = &self.vertices[triangle.vertices[1] as usize];
        let c = &self.vertices[triangle.vertices[2] as usize];

        let neighbour = &self.triangles[neighbour_idx as usize];
        let d = &self.vertices[neighbour.vertices[triangle.index_in_neighbours[2] as usize] as usize];

        let test = in_circle_2d(a.x_scaled, a.y_scaled, b.x_scaled, b.y_scaled, c.x_scaled, c.y_scaled, d.x_scaled, d.y_scaled);

        if test < 0. {
            self.flip_triangles(triangle_idx, neighbour_idx);
        }
    }

    fn flip_triangles(&mut self, triangle_idx: i32, neighbour_idx: i32) {
        let triangle = &self.triangles[triangle_idx as usize];
        let neighbour = &self.triangles[neighbour_idx as usize];

        // read necessary info
        let (ai, bi, ci) = (triangle.vertices[0], triangle.vertices[1], triangle.vertices[2]);
        let (na, nb) = (triangle.neighbours[0], triangle.neighbours[1]);
        let idx_in_na = triangle.index_in_neighbours[0];
        let idx_in_nb = triangle.index_in_neighbours[1];

        let idx_in_neighbour = triangle.index_in_neighbours[2] as usize;
        let fi = neighbour.vertices[idx_in_neighbour];
        let nd = neighbour.neighbours[(idx_in_neighbour + 1) % 3];
        let ne = neighbour.neighbours[(idx_in_neighbour + 2) % 3];
        let idx_in_nd = neighbour.index_in_neighbours[(idx_in_neighbour + 1) % 3];
        let idx_in_ne = neighbour.index_in_neighbours[(idx_in_neighbour + 2) % 3];

        // create new triangles
        let t0 = self.new_triangle_at(ai, fi, ci, triangle_idx);
        let t1 = self.new_triangle_at(fi, bi, ci, neighbour_idx);

        // fix neighbours
        self.triangles[t0 as usize].update_neighbours(t1, nb, nd, 1, idx_in_nb, idx_in_nd);
        self.triangles[t1 as usize].update_neighbours(na, t0, ne, idx_in_na, 0, idx_in_ne);
        self.triangles[na as usize].update_neighbour(t1, 0, idx_in_na);
        self.triangles[nd as usize].update_neighbour(t0, 2, idx_in_nd);
        // here the neighbours do not change, but idx_in_neighbour might change!
        self.triangles[nb as usize].update_neighbour(t0, 1, idx_in_nb);
        self.triangles[ne as usize].update_neighbour(t1, 2, idx_in_ne);

        self.triangles_to_check.push_back(t0);
        self.triangles_to_check.push_back(t1);

        // update current_triangle
        self.current_triangle_idx = t1;
    }

    fn consistency_check(&self) {
        if ! cfg!(debug_assertions) {
            return
        }
        // for each triangle check neighbour symmetry and check delaunay criterion
        // for each vertex check triangle information
        for (i, triangle) in self.triangles[3..].iter().enumerate() {
            let triangle_idx = (i + 3) as i32;
            for (j, &ngbr) in triangle.neighbours.iter().enumerate() {
                // check neighbouring relation symmetry
                let idx_in_ngbr = triangle.index_in_neighbours[j] as usize;
                assert_eq!(triangle_idx,
                           self.triangles[ngbr as usize].neighbours[idx_in_ngbr],
                           "Testing neighbour symmetry");

                // check overlapping vertices
                let neighbour = &self.triangles[ngbr as usize];
                let mut n_overlap = 0;
                for i in triangle.vertices.iter() {
                    if neighbour.vertices.contains(i) {
                        n_overlap += 1;
                    }
                }
                assert_eq!(n_overlap, 2, "Neighbours should have exactly 2 vertices in common!");

                // check Delaunay criterion
                // skip dummy triangles
                if ngbr < 3 {
                    continue;
                }
                let (a, b, c) = (&self.vertices[triangle.vertices[0] as usize],
                                 &self.vertices[triangle.vertices[1] as usize],
                                 &self.vertices[triangle.vertices[2] as usize]);
                let d = &self.vertices[neighbour.vertices[idx_in_ngbr] as usize];
                let in_circle = in_circle_2d(a.x_scaled, a.y_scaled, b.x_scaled, b.y_scaled, c.x_scaled, c.y_scaled, d.x_scaled, d.y_scaled);
                assert!(in_circle >= -1e-10, "in_circle test gave {:?} for vertices: \n a: {:?}\n b: {:?}\n c: {:?}\n d: {:?}\n", in_circle, a, b, c, d)
            }
        }
        for (i, vertex) in self.vertices.iter().enumerate() {
            assert_eq!(i as i32,
                       self.triangles[vertex.triangle as usize].vertices[vertex.index_in_triangle as usize],
                       "Testing vertex-triangle correspondence");
        }
        // println!("Consistency checks passed!");
    }

    pub fn to_str(&self) -> String {
        let mut result = String::from("# Vertices #\n");
        for (i, v) in self.vertices.iter().enumerate() {
            result += &format!("{}\t({}, {})\n", i, v.x, v.y);
        }

        result += "\n# Triangles #\n";
        for (i, triangle) in self.triangles[3..].iter().enumerate() {
            result += &format!("{}\t({}, {}, {})\n", i, triangle.vertices[0], triangle.vertices[1], triangle.vertices[2]);
        }
        result
    }

    pub fn to_file(&self, filename: &str) {
        fs::write(filename, self.to_str()).expect("Unable to write to file!");
    }
}