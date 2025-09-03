
file_path = "/Users/shengge/projetcs/robotics/kidy/urdf/AlphaBot_1S_B_RM75_description/urdf/scene.xml"

import time

import mujoco
import socket
import struct
import mujoco.viewer

m = mujoco.MjModel.from_xml_path(file_path)
d = mujoco.MjData(m)

cj1_id = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_ACTUATOR, "RM75_6F_r_joint1")
cj2_id = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_ACTUATOR, "RM75_6F_r_joint2")
cj3_id = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_ACTUATOR, "RM75_6F_r_joint3")
cj4_id = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_ACTUATOR, "RM75_6F_r_joint4")
cj5_id = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_ACTUATOR, "RM75_6F_r_joint5")
cj6_id = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_ACTUATOR, "RM75_6F_r_joint6")
cj7_id = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_ACTUATOR, "RM75_6F_r_joint7")

j1_id = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_JOINT, "RM75_6F_r_joint1")
j2_id = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_JOINT, "RM75_6F_r_joint2")
j3_id = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_JOINT, "RM75_6F_r_joint3")
j4_id = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_JOINT, "RM75_6F_r_joint4")
j5_id = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_JOINT, "RM75_6F_r_joint5")
j6_id = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_JOINT, "RM75_6F_r_joint6")
j7_id = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_JOINT, "RM75_6F_r_joint7")

# Example usage
HOST = "127.0.0.1"
PORT = 8080
BUFFER_SIZE = 1024

with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
    s.connect((HOST, PORT))
    with mujoco.viewer.launch_passive(m, d) as viewer:
    # Close the viewer automatically after 30 wall-seconds.
        while viewer.is_running():
            step_start = time.time()

            # mj_step can be replaced with code that also evaluates
            # a policy and applies a control signal before stepping the physics.
            mujoco.mj_step(m, d)

            # Pack the joint positions as 64-bit floats and send them to the server
            joint_positions = [d.qpos[j1_id], d.qpos[j2_id], d.qpos[j3_id], d.qpos[j4_id], d.qpos[j5_id], d.qpos[j6_id], d.qpos[j7_id]]
            data = struct.pack('7d', *joint_positions)  # 'd' for 64-bit float
            s.sendall(data)

            response = s.recv(BUFFER_SIZE)
            torques = struct.unpack('8d', response)
            # Apply the received torques to the actuators
            d.ctrl[cj1_id] = torques[1]
            d.ctrl[cj2_id] = torques[2]
            d.ctrl[cj3_id] = torques[3]
            d.ctrl[cj4_id] = torques[4]
            d.ctrl[cj5_id] = torques[5]
            d.ctrl[cj6_id] = torques[6]
            d.ctrl[cj7_id] = torques[7]
            # Example modification of a viewer option: toggle contact points every two seconds.
            with viewer.lock():
                viewer.opt.flags[mujoco.mjtVisFlag.mjVIS_CONTACTPOINT] = int(d.time % 2)

            # Pick up changes to the physics state, apply perturbations, update options from GUI.
            viewer.sync()

            # Rudimentary time keeping, will drift relative to wall clock.
            time_until_next_step = m.opt.timestep - (time.time() - step_start)
            if time_until_next_step > 0:
                time.sleep(time_until_next_step)