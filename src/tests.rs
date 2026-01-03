

fn main()
{
    use blanko_quaternions::quaternion::*;

    let qs = Quaternion::scaled_rotor(Angle::degrees(45.0), 0.0, 1.0, 1.0, 2.0);
    let qn = Quaternion::rotor(Angle::degrees(45.0), 0.0, 1.0, 1.0);

    println!("{:?}",  qs.transform_vector_scaled(&[1.0,0.0,0.0]) );
}
