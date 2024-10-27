pub fn assert_near<T: PartialOrd + std::ops::Sub<Output = T> + Copy + std::fmt::Display>(
    lhs: T,
    rhs: T,
    epsilon: T,
) -> () {
    let diff = if lhs > rhs { lhs - rhs } else { rhs - lhs };
    if diff > epsilon {
        println!("ERROR: assert_near({lhs}, {rhs}, {epsilon} failed.");
        println!("       diff = {diff}");
        assert!(false)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn valid_assert() {
        assert_near(1.1, 1., 0.10001);
    }

    #[test]
    #[should_panic]
    fn invalid_assert_throws() {
        assert_near(1.1, 1., 0.099);
    }
}
