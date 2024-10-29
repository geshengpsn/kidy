use crate::kinematics::KiChain;

struct Body {}

struct MultiBody {}

struct MultiBodyInfo {}

impl MultiBody {
    pub fn from_urdf() -> Self {
        todo!()
    }

    pub fn info() -> MultiBodyInfo {
        todo!()
    }

    pub fn get_chian<const N: usize>() -> KiChain<N> {
        todo!()
    }
}
