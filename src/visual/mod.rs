use std::{fs::File, io::BufReader};

use crate::model::utils::pose_to_se3;
use crate::Model;
use nalgebra::Matrix4;
use rerun::{
    components::Translation3D, datatypes::UVec3D, Arrows3D, Asset3D, Color, Mat3x3, Mesh3D,
    Position3D, Scale3D, Transform3D, TriangleIndices, Vec3D,
};
use urdf_rs::Vec3;

pub struct RerunVisualizer<'a> {
    pub rerun: rerun::RecordingStream,
    package_dir: String,
    model: &'a Model,
}

impl<'a> RerunVisualizer<'a> {
    pub fn new(id: &str, package_dir: &str, model: &'a Model) -> Self {
        let rerun = rerun::RecordingStreamBuilder::new(id).spawn().unwrap();
        rerun.set_time("timeline", std::time::SystemTime::now());
        let mut viewer = Self {
            rerun,
            model,
            package_dir: package_dir.to_string(),
        };
        for (i, link) in viewer.model.links.iter().enumerate() {
            viewer.add_link_visaulize_mesh(i, &link.global_zero_pose);
        }
        viewer
    }

    pub fn set_joints(&mut self, joints: &[f64]) {
        self.rerun
            .set_time("timeline", std::time::SystemTime::now());
        let links = self.model.fk(joints);
        for (i, link_pose) in links.iter().enumerate() {
            self.update_link_pose(link_pose, &self.model.links[i].name);
        }
    }

    fn update_link_pose(&mut self, pose: &Matrix4<f64>, link_name: &str) {
        let rot = pose.fixed_view::<3, 3>(0, 0);
        let trans = pose.fixed_view::<3, 1>(0, 3);
        // let rot_array = rot.as_array();
        let f32_array = [
            rot[(0, 0)] as f32,
            rot[(1, 0)] as f32,
            rot[(2, 0)] as f32,
            rot[(0, 1)] as f32,
            rot[(1, 1)] as f32,
            rot[(2, 1)] as f32,
            rot[(0, 2)] as f32,
            rot[(1, 2)] as f32,
            rot[(2, 2)] as f32,
        ];
        self.rerun
            .log(
                "/".to_string() + link_name,
                &Transform3D::from_translation_mat3x3(
                    Translation3D::new(trans[0] as f32, trans[1] as f32, trans[2] as f32),
                    Mat3x3::from(f32_array),
                ),
            )
            .unwrap();
    }

    fn add_link_visaulize_mesh(&mut self, index: usize, global_pose: &Matrix4<f64>) {
        let link = self.model.links.get(index).unwrap();
        let link_name = link.name.clone();
        self.update_link_pose(global_pose, &link_name);
        for visual in link.visual.iter() {
            let pose = pose_to_se3(&visual.origin);
            let rot = pose.fixed_view::<3, 3>(0, 0);
            let trans = pose.fixed_view::<3, 1>(0, 3);
            let f32_array = [
                rot[(0, 0)] as f32,
                rot[(1, 0)] as f32,
                rot[(2, 0)] as f32,
                rot[(0, 1)] as f32,
                rot[(1, 1)] as f32,
                rot[(2, 1)] as f32,
                rot[(0, 2)] as f32,
                rot[(1, 2)] as f32,
                rot[(2, 2)] as f32,
            ];
            if let urdf_rs::Geometry::Mesh { filename, scale } = &visual.geometry {
                let scale = scale.unwrap_or(Vec3([1.0, 1.0, 1.0]));
                let name = visual.name.clone().unwrap_or("visual".into());
                let entity_path = link_name.clone() + "/" + &name;

                // mesh

                let filename = if filename.starts_with("package://") {
                    filename.replace("package://", &self.package_dir)
                } else {
                    filename.clone()
                };
                println!("filename: {:?}", filename);

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
                    self.rerun
                        .log(
                            "/".to_string() + &entity_path,
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
                    self.rerun
                        .log(
                            "/".to_string() + &entity_path,
                            &Asset3D::from_file_path(filename).unwrap(),
                        )
                        .unwrap();
                }

                // transfrom
                self.rerun
                    .log(
                        "/".to_string() + &entity_path,
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

        if let Some(j) = link.joint.urdf_joint.as_ref() {
            let name = j.name.clone();
            let entity_path = link_name.clone() + "/" + &name;
            self.rerun
                .log(
                    "/".to_string() + &entity_path,
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
