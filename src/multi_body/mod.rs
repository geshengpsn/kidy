use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::BufReader,
    path::Path,
};

use liealg::{Algebra, SE3, SO3};
use petgraph::visit::Bfs;
use rerun::{
    components::Translation3D, datatypes::UVec3D, Arrows3D, Asset3D, Color, Mat3x3, Mesh3D,
    Position3D, RecordingStream, Scale3D, Transform3D, TriangleIndices, Vec3D,
};
use urdf_rs::{read_file, read_from_string, Vec3};

pub use urdf_rs;

#[derive(Debug)]
pub struct Link {
    pub joint: Option<urdf_rs::Joint>,
    pub joint_screw: Option<liealg::se3<f64>>,
    pub joint_relative_pose: liealg::SE3<f64>,
    pub joint_value: f64,
    pub global_pose: liealg::SE3<f64>,
    pub urdf_link: urdf_rs::Link,
}

#[derive(Debug)]
pub struct MultiBodyGraph {
    graph: petgraph::graphmap::DiGraphMap<usize, ()>,
    pub link_map: HashMap<usize, Link>,
    pub root_index: usize,
    pub leafs_index: Vec<usize>,
    pub name: String,
    pub material: Vec<urdf_rs::Material>,
}

impl MultiBodyGraph {
    pub fn bfs(&self, start: usize) -> Vec<usize> {
        let bfs = petgraph::visit::Bfs::new(&self.graph, start);
        let iter = BfsIter {
            graph: &self.graph,
            bfs,
        };
        iter.collect()
    }

    pub fn get_link(&self, index: usize) -> Option<&Link> {
        self.link_map.get(&index)
    }

    pub fn get_mut_link(&mut self, index: usize) -> &mut Link {
        self.link_map.get_mut(&index).unwrap()
    }

    pub fn parent(&self, index: usize) -> Option<usize> {
        self.graph
            .neighbors_directed(index, petgraph::Direction::Incoming)
            .next()
    }

    pub fn children(&self, index: usize) -> Vec<usize> {
        self.graph
            .neighbors_directed(index, petgraph::Direction::Outgoing)
            .collect()
    }

    pub fn set_joint_value(&mut self, joint_values: HashMap<String, f64>) {
        for (name, value) in joint_values {
            let link = self
                .link_map
                .values_mut()
                .filter(|link| link.joint.is_some())
                .find(|link| link.joint.as_ref().unwrap().name == name)
                .unwrap();
            link.joint_value = value;
        }
    }

    pub fn move_joint(&mut self) {
        let bfs = self.bfs(self.root_index);

        for link in bfs {
            let parent = self.parent(link);
            if let Some(parent_index) = parent {
                let parent_global_pose = self.get_link(parent_index).unwrap().global_pose.clone();
                let joint_value = self.get_link(link).unwrap().joint_value;
                let relative_pose = relative_pose(
                    self.get_link(link).unwrap().joint.as_ref().unwrap(),
                    joint_value,
                );
                self.get_mut_link(link).global_pose = parent_global_pose * relative_pose;
            }
        }
    }

    // pub fn get_chain(&self, start: &str, end: &str) -> Vec<usize> {
    //     let mut chain = vec![end];
    //     self.link_map
    //         .iter()
    //         .find(|(_, link)| link.urdf_link.name == end)
    //         .map(|(index, _)| *index)
    //         .map(|index| {
    //             chain.extend(self.get_chain_recursive(index, start));
    //         });
    //     loop {
    //         let parent = self.parent(end);
    //         match parent {
    //             Some(p_index) => {
    //                 chain.push(p_index);
    //                 if p_index == start {
    //                     break;
    //                 }
    //             }
    //             None => break,
    //         }
    //     }
    //     chain.reverse();
    //     chain
    // }
}

pub struct BfsIter<'a> {
    graph: &'a petgraph::graphmap::DiGraphMap<usize, ()>,
    bfs: Bfs<usize, HashSet<usize>>,
}

impl<'a> Iterator for BfsIter<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        self.bfs.next(self.graph)
    }
}

pub fn parse_urdf_string(string: &str) -> Result<MultiBodyGraph, urdf_rs::UrdfError> {
    let robot = read_from_string(string)?;
    Ok(parse_robot(robot))
}

pub fn parse_urdf_file(path: impl AsRef<Path>) -> Result<MultiBodyGraph, urdf_rs::UrdfError> {
    let robot = read_file(path)?;
    Ok(parse_robot(robot))
}

fn parse_robot(robot: urdf_rs::Robot) -> MultiBodyGraph {
    let map: HashMap<Option<String>, (usize, Option<urdf_rs::Joint>, urdf_rs::Link)> = robot
        .links
        .into_iter()
        .enumerate()
        .map(|(index, link)| {
            let j = robot
                .joints
                .iter()
                .find(|joint| joint.child.link == link.name);
            (j.map(|j| j.parent.link.clone()), (index, j.cloned(), link))
        })
        .collect::<HashMap<_, _>>();
    let root = map
        .iter()
        .find(|(parent, _)| parent.is_none())
        .map(|(_, (index, _, _))| *index)
        .expect("no root link");

    let leafs = map
        .iter()
        .filter(|(_, (_, _, link))| {
            robot
                .joints
                .iter()
                .all(|joint| joint.parent.link != link.name.clone())
        })
        .map(|(_, (index, _, _))| *index)
        .collect::<Vec<_>>();
    let mut graph = petgraph::graphmap::DiGraphMap::new();

    map.iter().for_each(|(_, (self_index, _, _))| {
        graph.add_node(*self_index);
    });

    map.iter().for_each(|(_, (self_index, _, link))| {
        if let Some((child_index, _, _)) = map.get(&Some(link.name.clone())) {
            graph.add_edge(*self_index, *child_index, ());
        }
    });

    let link_map = map
        .into_iter()
        .map(|(_, (i, j, l))| {
            (
                i,
                Link {
                    joint_value: 0.,
                    joint_screw: j.as_ref().and_then(joint_screw),
                    joint_relative_pose: j
                        .as_ref()
                        .map(joint_relative_pose)
                        .unwrap_or(liealg::SE3::<f64>::identity()),
                    global_pose: liealg::SE3::<f64>::identity(),
                    urdf_link: l,
                    joint: j,
                },
            )
        })
        .collect::<HashMap<_, _>>();

    let mut multi_body = MultiBodyGraph {
        graph,
        link_map,
        name: robot.name,
        material: robot.materials,
        root_index: root,
        leafs_index: leafs,
    };

    multi_body.move_joint();
    multi_body
}

fn relative_pose(joint: &urdf_rs::Joint, joint_value: f64) -> liealg::SE3<f64> {
    let pose = joint_relative_pose(joint);
    let screw = joint_screw(joint);
    if let Some(screw) = screw {
        pose * (screw * joint_value).exp()
    } else {
        SE3::identity()
    }
}

fn joint_relative_pose(joint: &urdf_rs::Joint) -> liealg::SE3<f64> {
    origin_to_se3(&joint.origin)
}

pub fn origin_to_se3(origin: &urdf_rs::Pose) -> liealg::SE3<f64> {
    let rpy = origin.rpy.0;
    let xyz = origin.xyz.0;
    SE3::new(&SO3::from_euler_angles(rpy[0], rpy[1], rpy[2]), xyz)
}

fn joint_screw(joint: &urdf_rs::Joint) -> Option<liealg::se3<f64>> {
    match joint.joint_type {
        urdf_rs::JointType::Revolute => {
            let norm =
                (joint.axis.xyz[0].powi(2) + joint.axis.xyz[1].powi(2) + joint.axis.xyz[2].powi(2))
                    .sqrt();
            let norm_x = joint.axis.xyz[0] / norm;
            let norm_y = joint.axis.xyz[1] / norm;
            let norm_z = joint.axis.xyz[2] / norm;
            Some(liealg::se3::<f64>::new(
                [norm_x, norm_y, norm_z],
                [0., 0., 0.],
            ))
        }
        _ => None,
    }
}

pub struct URDFViewer {
    recording_stream: RecordingStream,
    pub robot: MultiBodyGraph,
}

impl URDFViewer {
    pub fn new(urdf_file: impl AsRef<Path>) -> Self {
        let robot = parse_urdf_file(urdf_file).unwrap();
        let rr = rerun::RecordingStreamBuilder::new(robot.name.clone())
            .spawn()
            .unwrap();
        Self {
            recording_stream: rr,
            robot,
        }
    }
}

impl URDFViewer {
    pub fn visaulize_robot(&mut self) {
        let root = self.robot.root_index;
        let bfs = self.robot.bfs(root).to_vec();
        bfs.iter().for_each(|index| self.visaulize_link(*index));
    }

    pub fn transform_links(&mut self) {
        let bfs = self.robot.bfs(self.robot.root_index);
        for index in bfs {
            let link = self.robot.get_link(index).unwrap();
            let global_pose = &link.global_pose;
            let link_name = link.urdf_link.name.clone();
            let (rot, trans) = global_pose.rot_trans();
            let rot_array = rot.as_array();
            let f32_array = [
                rot_array[0] as f32,
                rot_array[1] as f32,
                rot_array[2] as f32,
                rot_array[3] as f32,
                rot_array[4] as f32,
                rot_array[5] as f32,
                rot_array[6] as f32,
                rot_array[7] as f32,
                rot_array[8] as f32,
            ];
            self.recording_stream
                .log(
                    link_name.clone(),
                    &Transform3D::from_translation_mat3x3(
                        Translation3D::new(trans[0] as f32, trans[1] as f32, trans[2] as f32),
                        Mat3x3::from(f32_array),
                    ),
                )
                .unwrap();
        }
    }

    fn visaulize_link(&mut self, index: usize) {
        let link = self.robot.get_link(index).unwrap();
        let global_pose = &link.global_pose;
        let link_name = link.urdf_link.name.clone();
        let (rot, trans) = global_pose.rot_trans();
        let rot_array = rot.as_array();
        let f32_array = [
            rot_array[0] as f32,
            rot_array[1] as f32,
            rot_array[2] as f32,
            rot_array[3] as f32,
            rot_array[4] as f32,
            rot_array[5] as f32,
            rot_array[6] as f32,
            rot_array[7] as f32,
            rot_array[8] as f32,
        ];
        self.recording_stream
            .log(
                link_name.clone(),
                &Transform3D::from_translation_mat3x3(
                    Translation3D::new(trans[0] as f32, trans[1] as f32, trans[2] as f32),
                    Mat3x3::from(f32_array),
                ),
            )
            .unwrap();
        for visual in link.urdf_link.visual.iter() {
            let pose = origin_to_se3(&visual.origin);
            let (rot, trans) = pose.rot_trans();
            let rot_array = rot.as_array();
            let f32_array = [
                rot_array[0] as f32,
                rot_array[1] as f32,
                rot_array[2] as f32,
                rot_array[3] as f32,
                rot_array[4] as f32,
                rot_array[5] as f32,
                rot_array[6] as f32,
                rot_array[7] as f32,
                rot_array[8] as f32,
            ];
            if let urdf_rs::Geometry::Mesh { filename, scale } = &visual.geometry {
                let scale = scale.unwrap_or(Vec3([1.0, 1.0, 1.0]));
                let name = visual.name.clone().unwrap_or("visual".into());
                let entity_path = link_name.clone() + "/" + &name;

                // mesh
                println!("filename: {:?}", filename);
                let filename = if filename.starts_with("package://") {
                    filename.replace("package://", "")
                } else {
                    filename.clone()
                };

                if filename.ends_with(".obj") {
                    let input = BufReader::new(File::open(&filename).unwrap());
                    let obj = obj::raw::object::parse_obj(input).unwrap();
                    let vertex = obj.positions.iter().map(|v| Position3D::new(v.0, v.1, v.2));
                    // let normals = obj.positions.iter().map(|v| {
                    //     Position3D::new(v.0, v.1, v.2)
                    // });
                    let indices = obj.polygons.iter().map(|p| {
                        if let obj::raw::object::Polygon::P(i) = p {
                            TriangleIndices(UVec3D::new(i[0] as u32, i[1] as u32, i[2] as u32))
                        } else {
                            panic!("only support triangle mesh")
                        }
                    });
                    let color = visual
                        .material
                        .as_ref()
                        .map(|m| {
                            color_name::Color::val()
                                .by_string(m.name.clone())
                                .unwrap_or(color_name::colors::white)
                        })
                        .unwrap_or(color_name::colors::white);
                    self.recording_stream
                        .log(
                            &*entity_path,
                            &Mesh3D::new(vertex)
                                .with_triangle_indices(indices)
                                // .with_vertex_normals(vertex_normals)
                                .with_vertex_colors(vec![
                                    Color::from_rgb(
                                        color[0], color[1], color[2]
                                    );
                                    obj.positions.len() * 3
                                ]),
                        )
                        .unwrap();
                } else {
                    self.recording_stream
                        .log(&*entity_path, &Asset3D::from_file(filename).unwrap())
                        .unwrap();
                }

                // transfrom
                self.recording_stream
                    .log(
                        entity_path,
                        &Transform3D::from_translation_mat3x3(
                            Translation3D::new(trans[0] as f32, trans[1] as f32, trans[2] as f32),
                            Mat3x3::from(f32_array),
                        )
                        .with_scale(Scale3D(Vec3D([
                            scale.0[0] as f32,
                            scale.0[1] as f32,
                            scale.0[2] as f32,
                        ]))),
                    )
                    .unwrap();
            }
        }

        if let Some(j) = link.joint.as_ref() {
            let name = j.name.clone();
            let entity_path = link_name.clone() + "/" + &name;
            self.recording_stream
                .log(
                    &*entity_path,
                    &Arrows3D::from_vectors(vec![Vec3D::new(
                        j.axis.xyz.0[0] as f32,
                        j.axis.xyz.0[1] as f32,
                        j.axis.xyz.0[2] as f32,
                    )]),
                )
                .unwrap();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_urdf_file() {
        let robot = parse_urdf_file("urdf/kuka_allegro_description/kuka.urdf").unwrap();
        robot.bfs(robot.root_index).iter().for_each(|j| {
            println!("{:?}", j);
        });

        println!("root index: {:?}", robot.root_index);
        let children = robot.children(robot.root_index);
        println!("root children: {:?}", children);

        println!("leaf index: {:?}", robot.leafs_index);
        let parent = robot.parent(robot.leafs_index[0]);
        println!("leaf parent: {:?}", parent);
    }
}
