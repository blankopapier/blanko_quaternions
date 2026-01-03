

fn main()
{
    test_dual_quat_exp_brute_force();
}




/// Check what the power series of exp() returns for dual quats
fn test_dual_quat_exp_brute_force()
{
    use blanko_quaternions::dual_quaternion::*;

    let dq_with_scalar    = DualQuaternion { w: 1.0, i: 2.0, j: 3.0, k: 4.0, ie: 5.0, je: 6.0, ke: 7.0, we: 8.0 };
    let dq_without_scalar = DualQuaternion { w: 0.0, i: 2.0, j: 3.0, k: 4.0, ie: 5.0, je: 6.0, ke: 7.0, we: 0.0 };

    let exp = |dq: DualQuaternion| -> DualQuaternion {
        let mut result = DualQuaternion { w: 1.0, i: 0.0, j: 0.0, k: 0.0, ie: 0.0, je: 0.0, ke: 0.0, we: 0.0 };
        let mut nth_power = DualQuaternion { w: 1.0, i: 0.0, j: 0.0, k: 0.0, ie: 0.0, je: 0.0, ke: 0.0, we: 0.0 };
        let mut fac = 1;

        // 10 should be enough, even 5 should do the trick but just to be sure...
        for i in 1..11
        {
            nth_power = nth_power * dq;
            fac *= i;

            let s = nth_power * (1.0 / (fac as f32));

            result += s;
        }

        result
    };

    println!("{}", exp(dq_with_scalar));
    println!("{}", exp(dq_without_scalar));
}
