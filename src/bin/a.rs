use kidy::Model;
use std::io::{Read, Write};
use std::net::{TcpListener, TcpStream};
use std::sync::Arc;
use std::thread;

fn main() {
    let model = Model::from_urdf("./urdf/rm_75_6fb_description/urdf/RM75-6F.urdf").unwrap();
    // let g = SO3::<f64>::from_euler_angles(0., FRAC_PI_3, FRAC_PI_2).inv().act(&Point::new(0., 0., -9.8));
    // let g = [8.4870489570875, 0., -4.9];
    // model.gravity_forces(q, gravity)
    start_server(Arc::new(model));
}

fn handle_client(mut stream: TcpStream, model: Arc<Model>) {
    let mut buffer = [0; 1024];
    while let Ok(size) = stream.read(&mut buffer) {
        if size == 0 {
            break;
        }
        let mut q: Vec<f64> = buffer[..size]
            .chunks(8)
            .take(7)
            .filter_map(|chunk| {
                if chunk.len() == 8 {
                    Some(f64::from_le_bytes(chunk.try_into().unwrap()))
                } else {
                    None
                }
            })
            .collect();
        q.insert(0, 0.);
        let gravity = [8.4870489570875, 0., -4.9];
        let forces = model.gravity_forces(&q, &gravity);
        let forces = model.get_actuator_torque(&forces);
        let response: Vec<u8> = forces.iter().flat_map(|f| f.to_le_bytes()).collect();
        stream.write_all(&response).unwrap();
    }
}

fn start_server(model: Arc<Model>) {
    let listener = TcpListener::bind("127.0.0.1:8080").unwrap();
    println!("Server listening on 127.0.0.1:8080");

    for stream in listener.incoming() {
        match stream {
            Ok(stream) => {
                let m = model.clone();
                thread::spawn(move || handle_client(stream, m));
            }
            Err(e) => eprintln!("Connection failed: {}", e),
        }
    }
}
