use std::{
    collections::{HashMap, HashSet},
    path::Path,
};

use liealg::prelude::*;
use nalgebra::{Matrix3, Matrix6};
use petgraph::visit::Bfs;
use urdf_rs::read_file;

#[derive(Debug, Clone)]
struct Inertia {
    mass: f64,
    inertia: RotateInertia,
}

#[derive(Debug, Clone)]
struct RotateInertia {
    xx: f64,
    yy: f64,
    zz: f64,
    xy: f64,
    xz: f64,
    yz: f64,
}

impl Inertia {
    fn from_link(link: &urdf_rs::Link) -> Inertia {
        Inertia::new(
            link.inertial.mass.value,
            RotateInertia {
                xx: link.inertial.inertia.ixx,
                yy: link.inertial.inertia.iyy,
                zz: link.inertial.inertia.izz,
                xy: link.inertial.inertia.ixy,
                xz: link.inertial.inertia.ixz,
                yz: link.inertial.inertia.iyz,
            },
        )
    }

    fn new(mass: f64, inertia: RotateInertia) -> Inertia {
        Inertia { mass, inertia }
    }

    fn spatial_inertia(&self) -> Matrix6<f64> {
        let mut res = Matrix6::zeros();
        res.fixed_view_mut::<3, 3>(0, 0).copy_from(&Matrix3::new(
            self.inertia.xx,
            self.inertia.xy,
            self.inertia.xz,
            self.inertia.xy,
            self.inertia.yy,
            self.inertia.yz,
            self.inertia.xz,
            self.inertia.yz,
            self.inertia.zz,
        ));
        res.fixed_view_mut::<3, 3>(3, 3)
            .copy_from(&Matrix3::from_diagonal_element(self.mass));
        res
    }
}

#[derive(Debug, Clone)]
pub struct Link {
    pub space_screw: Option<liealg::se3<f64>>,
    pub local_screw: Option<liealg::se3<f64>>,
    pub local_zero_pose: liealg::SE3<f64>,
    pub global_zero_pose: liealg::SE3<f64>,
    pub local_spatial_inertial: Matrix6<f64>,

    pub joint: Option<urdf_rs::Joint>,
    pub urdf_link: urdf_rs::Link,
}

#[derive(Debug, Clone)]
pub struct MultiBody {
    // link index graph
    graph: petgraph::graphmap::DiGraphMap<usize, ()>,
    // map index -> link
    pub link_map: HashMap<usize, Link>,
    pub root_index: usize,
    pub leafs_index: Vec<usize>,
    pub name: String,
}

impl MultiBody {
    pub fn from_urdf(path: impl AsRef<Path>) -> Result<MultiBody, urdf_rs::UrdfError> {
        let robot = read_file(path)?;
        Ok(parse_robot(robot))
    }

    fn bfs(&self, start: usize) -> Vec<usize> {
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

    fn get_mut_link(&mut self, index: usize) -> &mut Link {
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

    pub fn get_kidy_chain<const N: usize>(&self, start: &str, end: &str) -> KidyChain<N> {
        let mut chain = vec![];
        let start = self
            .link_map
            .iter()
            .find(|(_, link)| link.urdf_link.name == start)
            .map(|(index, _)| *index)
            .unwrap();
        let end = self
            .link_map
            .iter()
            .find(|(_, link)| link.urdf_link.name == end)
            .map(|(index, _)| *index)
            .unwrap();

        let mut current = end;
        loop {
            chain.push(current);
            if let Some(p) = self.parent(current) {
                if p == start {
                    chain.push(p);
                    break;
                }
                current = p;
            } else {
                break;
            }
        }
        chain.reverse();
        let chain: Vec<_> = chain
            .iter()
            .map(|index| self.get_link(*index).unwrap().clone())
            .collect();
        let joints_screw: Vec<_> = chain
            .iter()
            .filter_map(|link| link.space_screw.clone())
            .collect();
        assert_eq!(joints_screw.len(), N);
        let zero_pose = chain[0].global_zero_pose.clone();
        let local_screw = chain
            .iter()
            .filter_map(|link| link.local_screw.clone())
            .collect();
        // assert_eq!(local_screw.len(), N);
        let mut local_zero_pose: Vec<_> = chain
            .iter()
            .filter(|l| match &l.joint {
                Some(j) => j.joint_type != urdf_rs::JointType::Fixed,
                None => true,
            })
            .map(|link| link.local_zero_pose.clone())
            .collect();
        local_zero_pose.push(chain.last().unwrap().local_zero_pose.clone());

        let local_spatial_inertial = chain
            .iter()
            .filter(|l| match &l.joint {
                Some(j) => j.joint_type != urdf_rs::JointType::Fixed,
                None => true,
            })
            .map(|link| link.local_spatial_inertial)
            .collect();

        KidyChain {
            joints_screw,
            zero_pose,
            local_screw,
            local_zero_pose,
            local_spatial_inertial,
        }
    }
}

struct BfsIter<'a> {
    graph: &'a petgraph::graphmap::DiGraphMap<usize, ()>,
    bfs: Bfs<usize, HashSet<usize>>,
}

impl<'a> Iterator for BfsIter<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        self.bfs.next(self.graph)
    }
}

fn parse_robot(robot: urdf_rs::Robot) -> MultiBody {
    // parent link name -> (link index, joint, link)
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

    // find root link index
    let root = map
        .iter()
        .find(|(parent, _)| parent.is_none())
        .map(|(_, (index, _, _))| *index)
        .expect("no root link");

    // find leaf link index
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
                    space_screw: None,
                    local_screw: None,
                    global_zero_pose: liealg::SE3::identity(),
                    local_zero_pose: liealg::SE3::identity(),
                    local_spatial_inertial: Inertia::from_link(&l).spatial_inertia(),
                    urdf_link: l,
                    joint: j,
                },
            )
        })
        .collect::<HashMap<_, _>>();

    let mut multi_body = MultiBody {
        graph,
        link_map,
        name: robot.name,
        root_index: root,
        leafs_index: leafs,
    };

    let bfs = multi_body.bfs(root);
    for link in bfs {
        let parent = multi_body.parent(link);
        if let Some(parent_index) = parent {
            let parent_global_pose = multi_body
                .get_link(parent_index)
                .unwrap()
                .global_zero_pose
                .clone();
            let relative_pose =
                joint_relative_pose(multi_body.get_link(link).unwrap().joint.as_ref().unwrap());
            multi_body.get_mut_link(link).local_zero_pose = relative_pose.clone();

            let global_pose = parent_global_pose * relative_pose;
            multi_body.get_mut_link(link).global_zero_pose = global_pose.clone();

            if let Some(screw) =
                joint_screw(multi_body.get_link(link).unwrap().joint.as_ref().unwrap())
            {
                multi_body.get_mut_link(link).local_screw = Some(screw.clone());
                let global_screw = global_pose.adjoint().act(&screw);
                multi_body.get_mut_link(link).space_screw = Some(global_screw);
            }
        }
    }
    multi_body
}

fn joint_screw(joint: &urdf_rs::Joint) -> Option<liealg::se3<f64>> {
    match joint.joint_type {
        urdf_rs::JointType::Revolute | urdf_rs::JointType::Continuous => {
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

fn joint_relative_pose(joint: &urdf_rs::Joint) -> liealg::SE3<f64> {
    origin_to_se3(&joint.origin)
}

fn origin_to_se3(origin: &urdf_rs::Pose) -> liealg::SE3<f64> {
    let rpy = origin.rpy.0;
    let xyz = origin.xyz.0;
    liealg::SE3::new(&liealg::SO3::from_euler_angles(rpy[0], rpy[1], rpy[2]), xyz)
}

#[derive(Debug)]
pub struct KidyChain<const N: usize> {
    // pub(crate) value: Vec<Link>,
    pub(crate) joints_screw: Vec<liealg::se3<f64>>,
    pub(crate) zero_pose: liealg::SE3<f64>,
    pub(crate) local_screw: Vec<liealg::se3<f64>>,

    pub(crate) local_zero_pose: Vec<liealg::SE3<f64>>,
    pub(crate) local_spatial_inertial: Vec<Matrix6<f64>>,
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn test_inertia() {
        let inertia = Inertia::new(
            1.0,
            RotateInertia {
                xx: 1.0,
                yy: 2.0,
                zz: 3.0,
                xy: 4.0,
                xz: 5.0,
                yz: 6.0,
            },
        );
        let spatial_inertia = inertia.spatial_inertia();
        println!("{}", spatial_inertia);
        assert_eq!(spatial_inertia[(0, 0)], 1.0);
        assert_eq!(spatial_inertia[(1, 1)], 2.0);
        assert_eq!(spatial_inertia[(2, 2)], 3.0);
        assert_eq!(spatial_inertia[(3, 3)], 1.0);
        assert_eq!(spatial_inertia[(4, 4)], 1.0);
        assert_eq!(spatial_inertia[(5, 5)], 1.0);
    }
    #[test]
    fn test_from_urdf() {
        let multi_body =
            MultiBody::from_urdf("urdf/franka_description/robots/franka_panda.urdf").unwrap();
        let chain = multi_body.get_kidy_chain::<7>("panda_link0", "panda_hand");

        // N
        // for p in chain.joints_screw {
        //     println!("{:.3}", p.vee());
        // }

        // N
        // for p in chain.local_screw {
        //     println!("{:.3}", p.vee());
        // }

        // N + 1
        for p in chain.local_zero_pose {
            println!("{:.3}", p);
        }

        // N + 1
        for p in chain.local_spatial_inertial {
            println!("{:.3}", p);
        }
        // println!("{:.3}", chain.zero_ee_pose());
    }
}
