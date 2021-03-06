use super::delaunay::DelaunayTriangulation;
use super::geometry::{Triangle, line_segment_intersection};
use super::NeighbourDirection;
use std::fs;
use std::fmt;
use std::iter::FromIterator;
use crate::vector::Vec2;
use std::fmt::Formatter;

#[derive(Copy, Clone)]
pub(super) struct IntersectionResult<'a> {
    pub(super) distance: f64,
    pub(super) intersection: Vec2<f64>,
    pub(super) exit_face_idx: i32,
    pub(super) exit_face: &'a VoronoiFace
}


#[derive(Debug, Clone)]
pub struct VoronoiError {
    pub message: String,
    pub line: u32,
    pub column: u32
}

impl fmt::Display for VoronoiError {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{} ({}:{})", self.message, self.line, self.column)
    }
}

impl std::error::Error for VoronoiError {}

/// A face (line) between two cells in a voronoi grid
#[derive(Debug, Default)]
pub(super) struct VoronoiFace {
    area: f64,
    start: Vec2<f64>,
    end: Vec2<f64>,
    midpoint: Vec2<f64>,
    pub(super) adjacent_cells: [i32; 2],
    pub(super) neighbour_direction: Option<NeighbourDirection>
}


/// A cell from a voronoi grid in 2D
#[derive(Debug)]
pub(super) struct VoronoiCell {
    vertices: Vec<i32>,
    faces: Vec<i32>,
    pub(super) centroid: Vec2<f64>,
    volume: f64
}

impl Default for VoronoiCell {
    fn default() -> VoronoiCell {
        VoronoiCell {
            vertices: Vec::new(),
            faces: Vec::new(),
            centroid: Vec2::new(0.0, 0.0),
            volume: 0.
        }
    }
}

#[derive(Default, Debug)]
pub struct VoronoiGrid {
    vertices: Vec<Vec2<f64>>,
    faces: Vec<VoronoiFace>,
    pub(super) cells: Vec<VoronoiCell>,
    is_periodic: bool,
    anchor: Vec2<f64>,
    sides: Vec2<f64>,
    n_cells: usize
}

impl VoronoiGrid {
    pub fn with_capacity(n_cells: usize, n_vertices: usize) -> VoronoiGrid {
        VoronoiGrid {
            vertices: Vec::with_capacity(n_vertices),
            faces: Vec::with_capacity(2 * n_vertices),
            cells: Vec::with_capacity(n_cells),
            ..VoronoiGrid::default()
        }
    }

    pub fn from_delaunay_triangulation(triangulation: &DelaunayTriangulation) -> VoronoiGrid {
        let mut grid = VoronoiGrid::with_capacity(triangulation.vertices.len() - 3,
                                                  triangulation.triangles.len() - 3);
        grid.n_cells = triangulation.n_vertices;
        grid.is_periodic = triangulation.is_periodic;
        assert!(grid.is_periodic, "Non periodic grids are not yet supported!");
        grid.anchor = triangulation.anchor + triangulation.sides;
        grid.sides = triangulation.sides;
        // for each triangle of triangulation add the circumcenter to vertices (skip dummy triangles)
        for triangle in triangulation.triangles[3..].iter() {
            grid.vertices.push(triangle.circumcenter(triangulation));
        }

        // For each vertex of triangulation (generator): loop around the Delaunay triangles containing
        // that vertex. The circumcenters of those triangles are the voronoi vertices of the cell
        // generated by the current generator. Update the area and centroid of the current cell as
        // you go.
        for i in 3..(triangulation.n_vertices + 3) {
            grid.add_cell_from_delaunay_generator(i, triangulation);
            // if triangulation.is_connected_to_non_dummy_non_ghost_vertex(i){
            //     grid.add_cell_from_delaunay_generator(i, triangulation);
            // }
        }
        grid
    }

    pub fn from_points(points_x: &Vec<f64>,
                       points_y: &Vec<f64>,
                       anchor: Vec2<f64>,
                       sides: Vec2<f64>,
                       make_periodic: bool) -> VoronoiGrid {
        let delaunay = DelaunayTriangulation::from_points(
            points_x,
            points_y,
            anchor,
            sides,
            make_periodic
        );
        VoronoiGrid::from_delaunay_triangulation(&delaunay)
    }

    fn add_cell_from_delaunay_generator(&mut self, generator_idx: usize, triangulation: &DelaunayTriangulation) {
        let generator = &triangulation.vertices[generator_idx];
        let current_voronoi_cell_idx = generator_idx as i32 - 3;
        let mut idx_in_current_triangle = generator.index_in_triangle;
        let mut current_cell = VoronoiCell::default();
        let generator_as_vertex2d = Vec2::new(generator.x, generator.y);

        for current_triangle_idx_in_d in triangulation.get_triangle_idx_around_vertex(generator_idx) {
            let current_triangle = &triangulation.triangles[current_triangle_idx_in_d];
            assert_eq!(current_voronoi_cell_idx + 3, current_triangle.vertices[idx_in_current_triangle as usize]);
            let next_triangle_idx_in_current_triangle = ((idx_in_current_triangle + 1) % 3) as usize;
            let next_triangle_idx_in_d = current_triangle.neighbours[next_triangle_idx_in_current_triangle];

            let current_voronoi_vertex_idx = current_triangle_idx_in_d as i32 - 3;
            let next_voronoi_vertex_idx: i32 = next_triangle_idx_in_d - 3;
            current_cell.vertices.push(current_voronoi_vertex_idx);

            if generator_idx < triangulation.n_vertices + 3 {
                let current_wedge = Triangle::new(
                    generator_as_vertex2d,
                    self.vertices[current_voronoi_vertex_idx as usize],
                    self.vertices[next_voronoi_vertex_idx as usize]
                );
                let current_wedge_area = current_wedge.area();
                // create faces between cells
                let neighbouring_generator_idx_in_d = current_triangle.vertices[((idx_in_current_triangle + 2) % 3) as usize];
                let neighbour_direction;
                let neighbouring_voronoi_cell_idx;
                if neighbouring_generator_idx_in_d < triangulation.n_vertices as i32 + 3 {
                    neighbour_direction = None;
                    neighbouring_voronoi_cell_idx = neighbouring_generator_idx_in_d - 3;
                } else {
                    let info = triangulation.neighbour_info[neighbouring_generator_idx_in_d as usize - triangulation.n_vertices - 3];
                    neighbour_direction = Some(info.1);
                    neighbouring_voronoi_cell_idx = info.0 as i32 - 3;
                };
                current_cell.faces.push(
                    self.create_or_get_face(
                        current_voronoi_vertex_idx,
                        next_voronoi_vertex_idx,
                        current_voronoi_cell_idx,
                        neighbouring_voronoi_cell_idx,
                        neighbour_direction
                    )
                );

                // Update area and area weighted centroid sum of cell
                current_cell.volume += current_wedge_area;
                current_cell.centroid += current_wedge.centroid() * current_wedge_area;
            }
            // update idx_in_current_triangle
            let current_triangle_idx_in_next_triangle = current_triangle.index_in_neighbours[next_triangle_idx_in_current_triangle];
            idx_in_current_triangle = (current_triangle_idx_in_next_triangle + 1) % 3;
        }
        // Divide area weighted sum of centroids by total area of cell -> centroid of cell
        current_cell.centroid /= current_cell.volume;
        self.cells.push(current_cell);
    }

    fn create_or_get_face(&mut self, vertex_from_idx: i32, vertex_to_idx: i32, cell_in_idx: i32, cell_out_idx: i32, neighbour_direction: Option<NeighbourDirection>) -> i32 {
        // assert_ne!(cell_in_idx, cell_out_idx, "Trying to add face between a cell and itself!");
        assert_ne!(vertex_from_idx, vertex_to_idx, "Trying to add a face from a vertex to itself!");
        let face_idx: i32;
        if cell_out_idx < cell_in_idx && cell_out_idx > 2 && neighbour_direction.is_none() {
            let cell_out = &self.cells[cell_out_idx as usize];
            let face_idx_in_cell_out = cell_out.vertices.iter().position(|&v_idx| v_idx == vertex_to_idx).unwrap();
            face_idx = cell_out.faces[face_idx_in_cell_out];
            assert_eq!(&self.faces[face_idx as usize].adjacent_cells[1], &cell_in_idx,
                       "Face {:?} has wrong adjacent cells!", self.faces[face_idx as usize]);
        } else {
            // create new face and return index
            face_idx = self.faces.len() as i32;
            self.faces.push(
                VoronoiFace {
                    area: (self.vertices[vertex_to_idx as usize] - self.vertices[vertex_from_idx as usize]).norm(),
                    start: self.vertices[vertex_from_idx as usize],
                    end: self.vertices[vertex_to_idx as usize],
                    midpoint: (self.vertices[vertex_to_idx as usize] + self.vertices[vertex_from_idx as usize]) / 2.,
                    adjacent_cells: [cell_in_idx, cell_out_idx],
                    neighbour_direction
                }
            );
        }
        face_idx
    }

    pub(super) fn raytrace_cell(&self, cell_idx: usize, start_x: f64, start_y: f64, cos_theta: f64, sin_theta: f64, prev_face_idx: i32) -> Result<IntersectionResult, VoronoiError> {
        let cell = &self.cells[cell_idx];
        let end_x = start_x + sin_theta * self.sides.x() * 2.;
        let end_y = start_y - cos_theta * self.sides.y() * 2.;

        let mut face;

        let mut possible_results = vec![];
        for &face_idx in cell.faces.iter() {
            if face_idx == prev_face_idx {
                continue;
            }
            face = &self.faces[face_idx as usize];
            let f_start_x = face.start.x();
            let f_start_y = face.start.y();
            let f_end_x = face.end.x();
            let f_end_y = face.end.y();
            if let Some((intersection_x, intersection_y)) = line_segment_intersection(start_x, start_y, end_x, end_y, f_start_x, f_start_y, f_end_x, f_end_y) {
                let dist = (Vec2::new(start_x, start_y) - Vec2::new(intersection_x, intersection_y)).norm();
                if dist > 1e-13 {
                    return Ok(IntersectionResult{
                        distance: dist,
                        intersection: Vec2::new(intersection_x, intersection_y),
                        exit_face_idx: face_idx,
                        exit_face: face
                    });
                } else {
                    possible_results.push(IntersectionResult{
                        distance: dist,
                        intersection: Vec2::new(intersection_x, intersection_y),
                        exit_face_idx: face_idx,
                        exit_face: face
                    })
                }
            }
        }
        if possible_results.len() == 0 {
            Err(VoronoiError{
                message: "No valid intersected face found!".to_string(), line: line!(), column: column!()
            })
        } else {
            let mut idx_max = 0;
            let mut cur_max = 0.;
            for (i, res) in possible_results.iter().enumerate() {
                if res.distance > cur_max {
                    idx_max = i;
                    cur_max = res.distance;
                }
            }
            Ok(possible_results[idx_max])
        }

    }

    pub fn lloyd_relax(&self, move_threshold: f64, max_iter: usize) -> VoronoiGrid {
        let mut generators_x = Vec::from_iter(self.cells[..self.n_cells].iter().map(|c| c.centroid.x()));
        let mut generators_y = Vec::from_iter(self.cells[..self.n_cells].iter().map(|c| c.centroid.y()));

        let mut d = DelaunayTriangulation::from_points(
            &generators_x,
            &generators_y,
            self.anchor,
            self.sides,
            self.is_periodic
        );
        let mut v: VoronoiGrid = VoronoiGrid::from_delaunay_triangulation(&d);

        let mut displacement_threshold_satisfied = false;
        let mut previous_max_displacement = f64::INFINITY;
        let mut iter: usize = 0;
        while iter < max_iter && !displacement_threshold_satisfied {
            let mut max_displacement = 0.;
            for (i, voronoi_cell) in v.cells[..v.n_cells].iter().enumerate() {
                let centroid = voronoi_cell.centroid;
                let generator = Vec2::new(generators_x[i], generators_y[i]);
                let displacement = (generator - centroid).norm();
                if displacement > max_displacement {
                    max_displacement = displacement
                }
            }
            if max_displacement > previous_max_displacement { break; }
            previous_max_displacement = max_displacement;
            displacement_threshold_satisfied = max_displacement < move_threshold;
            iter += 1;
            generators_x = Vec::from_iter(v.cells[..v.n_cells].iter().map(|c| c.centroid.x()));
            generators_y = Vec::from_iter(v.cells[..v.n_cells].iter().map(|c| c.centroid.y()));
            d = DelaunayTriangulation::from_points(
                &generators_x,
                &generators_y,
                self.anchor,
                self.sides,
                self.is_periodic
            );
            v = VoronoiGrid::from_delaunay_triangulation(&d);
            println!("Relaxation iter: {}, maximum displacement was: {}", iter, max_displacement);
        }
        v
    }

    pub fn to_str(&self) -> String {
        let mut result = String::from("# Vertices #\n");
        for (i, v) in self.vertices.iter().enumerate() {
            result += &format!("{}\t({}, {})\n", i, v.x(), v.y());
        }

        result += "\n# Cells #\n";
        for (i, cell) in self.cells.iter().enumerate() {
            result += &format!("{}\t(", i);
            for (i, vertex_idx) in cell.vertices.iter().enumerate(){
                result += &format!("{}", vertex_idx);
                if i < cell.vertices.len() - 1 {
                    result += ", "
                }
            }
            result += ")\n";
        }

        result += "\n# Centroids #\n";
        for (i, cell) in self.cells.iter().enumerate() {
            result += &format!("{}\t({}, {})\n", i, cell.centroid.x(), cell.centroid.y());
        }
        result
    }

    pub fn to_file(&self, filename: &str) {
        fs::write(filename, self.to_str()).expect("Unable to write to file!");
    }
}