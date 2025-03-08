use crate::spacepoint::SpacePoint;

#[derive(Debug, PartialEq)]
struct Doublet<'a>(pub &'a SpacePoint, pub &'a SpacePoint);

fn find_doublets<'a>(spacepoints: &'a Vec<SpacePoint>) -> Vec<Doublet<'a>> {
    if let Some(result) = spacepoints
        .iter()
        .map(|sp| {
            spacepoints
                .iter()
                .filter(move |sp2| sp != *sp2 && sp.is_compatible(sp2))
                .map(move |sp2| Doublet(sp, sp2))
                .collect::<Vec<Doublet>>()
        })
        .reduce(|mut v1, mut v2| {
            v1.append(&mut v2);
            v1
        })
    {
        return result;
    }
    return Vec::new();
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::matrix::Vector3;

    fn new_sp(x: f64, y: f64, z: f64) -> SpacePoint {
        Vector3::new(x, y, z).into()
    }

    #[test]
    fn find_doublets_() {
        let spacepoints = vec![
            new_sp(0., 0., 0.),
            new_sp(1., 2., 3.),
            new_sp(1., 2., 3.01),
            new_sp(1., 2.05, 3.),
            new_sp(0.01, 0., 0.01),
            new_sp(0., 0.01, 0.),
        ];

        let doublets = find_doublets(&spacepoints);

        let mut expected_doublets = Vec::new();
        for sp in &spacepoints {
            for sp2 in &spacepoints {
                if sp2 != sp && sp.is_compatible(sp2) {
                    expected_doublets.push(Doublet(&sp, &sp2));
                }
            }
        }

        assert!(
            expected_doublets.len() >= 5,
            "Use more significant test data."
        );

        assert_eq!(doublets.len(), expected_doublets.len());
        assert_eq!(doublets, expected_doublets);
    }
}
