

fn main()
{
    use blanko_quaternions::dual_quaternion::*;

    let dq = DualQuaternion::screw(
        &DualQuaternion::line(&[0.0,1.0,0.0], &[1.0,0.0,0.0]),
        Angle::degrees(90.0),
        1.0
    );

    let r1 = DualQuaternion::rotor(Angle::degrees(45.0), &[0.0,0.0,1.0]);
    let t1 = DualQuaternion::translator(&[0.0,1.0,0.0]);


    let p = (r1*dq).transform_point(&[0.0,0.0,0.0]);
    let l = DualQuaternion::line(&[0.0,0.0,0.0], &[1.0,0.0,0.0]);

    println!("{:?}", dq * l * dq.conj() );
}
