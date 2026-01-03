

fn main()
{
    //test_slerp_quat();
    //test_dual_quat_exp_brute_force();
    //test_dual_quat_log();
    //test_dual_quat_motor();
    test_dual_quat_sclerp();
}




fn test_slerp_quat()
{
    use blanko_quaternions::quaternion::*;

    let q1 = Quaternion::scaled_rotor(Angle::degrees(90.0), &[1.0, 0.0, 0.0], 0.5);
    let q2 = Quaternion::scaled_rotor(Angle::degrees(90.0), &[0.0, 0.0, 1.0], 2.0);

    println!("{:?}", (q1.slerp(q2, 1.0)).transform_vector_scaled(&[1.0,0.0,0.0]));
}

fn test_dual_quat_sclerp()
{
    use blanko_quaternions::dual_quaternion::*;

    let line1 = DualQuaternion::line(&[0.0,1.0,0.0], &[1.0,0.0,0.0]);
    let line2 = DualQuaternion::line(&[1.0,0.0,0.0], &[0.0,0.0,1.0]);

    let screw1 = DualQuaternion::screw(&line1, Angle::degrees(90.0), 2.0);
    let screw2 = DualQuaternion::screw(&line2, Angle::degrees(90.0), 3.0);

    let screw = screw1.sclerp(&screw2, 0.0);

    println!("{:?}", screw.transform_point(&[0.0,0.0,0.0]));
}



fn test_dual_quat_motor()
{
    use blanko_quaternions::dual_quaternion::*;

    let line1 = DualQuaternion::line(&[0.0,1.0,0.0], &[1.0,0.0,0.0]);
    let line2 = DualQuaternion::line(&[1.0,0.0,0.0], &[0.0,0.0,1.0]);

    let screw1 = DualQuaternion::screw(&line1, Angle::degrees(90.0), 2.0);
    let screw2 = DualQuaternion::screw(&line2, Angle::degrees(90.0), 3.0);

    println!("{:?}", (screw2*screw1).transform_point(&[0.0,0.0,0.0]));
}


/// Check whether or not the logarithm behaves correctly
fn test_dual_quat_log()
{
    use blanko_quaternions::dual_quaternion::*;

    let dq = DualQuaternion { w: 0.0, i: 2.0, j: 3.0, k: 4.0, ie: 5.0, je: 6.0, ke: 7.0, we: 0.0 };

    let dq_exp = dq.exp();
    let dq_log = dq_exp.log();

    println!("{:?}", dq_log.exp() );
    println!("{:?}", dq_exp );
}


/// Check what the power series of exp() returns for dual quats
fn test_dual_quat_exp_brute_force()
{
    use blanko_quaternions::dual_quaternion::*;

    let dq = DualQuaternion { w: 0.0, i: 2.0, j: 3.0, k: 4.0, ie: 5.0, je: 6.0, ke: 7.0, we: 0.0 };

    // Approximate exp(dq)
    let exp = |dq: DualQuaternion, max_i: usize| -> DualQuaternion {
        let mut result = DualQuaternion { w: 1.0, i: 0.0, j: 0.0, k: 0.0, ie: 0.0, je: 0.0, ke: 0.0, we: 0.0 };
        let mut nth_power = DualQuaternion { w: 1.0, i: 0.0, j: 0.0, k: 0.0, ie: 0.0, je: 0.0, ke: 0.0, we: 0.0 };
        let mut fac: usize = 1;

        for i in 1..max_i
        {
            nth_power = nth_power * dq;
            fac *= i;

            let s = nth_power * (1.0 / (fac as f32));

            result += s;
        }

        result
    };

    // Check if the approximation converges

    for i in 2..21
    {
        let r_margin = ( exp(dq, i) - dq.exp() ).norm();
        let i_margin = ( exp(dq, i) - dq.exp() ).inorm();

        println!("r: {}\ni: {}\n", r_margin, i_margin);
    }
}
