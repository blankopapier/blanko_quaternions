use blanko_quaternions::{direction, point};


fn main()
{
    use blanko_quaternions::dual_quaternion::*;

    let dq = DualQuaternion::from_line(
        &Point { x: 0.0, y: 1.0, z: 0.0 },
        &Direction { x: 1.0, y: 0.0, z: 0.0 },
        Angle::from_deg(90.0),
        1.0
    );

    let r1 = DualQuaternion::from_angle_axis(Angle::from_deg(45.0), &Direction { x: 0.0, y: 0.0, z: 1.0 });
    let t1 = DualQuaternion::from_translation(&Direction { x: 1.0, y: 0.0, z: 0.0 });


    let p = (r1*dq).transform_point(&Point { x: 0.0, y: 0.0, z: 0.0 });

    println!("{:?}", p);
}

fn test()
{
    use blanko_quaternions::point::*;

    point!(1,2,3);
}
