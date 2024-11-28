# kidy is a library for kinematics and dynamics of multi-body.

kidy = kinematics + dynamics

## get started
```rust
use kidy::*;

fn main() {
    let robot = MultiBody::from_urdf("path/to/urdf").unwrap();
    let chain = multi_body.get_kidy_chain::<7>("base_link0", "robot_hand");
}

```