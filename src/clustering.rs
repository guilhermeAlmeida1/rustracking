use crate::detector_module::*;

pub struct Hit {
    pub module_id: u64,
    pub pos: PixelPosition,
}

impl Hit {
    pub fn new(module_id: u64, pos: PixelPosition) -> Self {
        Hit { module_id, pos }
    }

    pub fn adjacent(&self, other: &Hit) -> bool {
        if self.module_id != other.module_id {
            return false;
        }

        let mut a = [self.pos.0, other.pos.0];
        let mut b = [self.pos.1, other.pos.1];
        a.sort();
        b.sort();

        if a[1] - a[0] > 1 || b[1] - b[0] > 1 {
            return false;
        }
        true
    }
}

pub fn ccl(hits: Vec<Hit>) -> Vec<usize> {
    let mut indices = vec![0; hits.len()];
    for i in 0..indices.len() {
        indices[i] = i;
    }
    for i in 0..indices.len() {
        for j in (0..i).rev() {
            if hits[i].adjacent(&hits[j]) {
                println!("Found adjacent: {i} {j}");
                indices[i] = j;
                break;
            }
        }
    }
    for i in 0..indices.len() {
        if indices[i] != i {
            let mut j = indices[i];
            while indices[j] != j {
                j = indices[j];
            }
            indices[i] = j;
        }
    }

    indices
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn adjacent() {
        let centre = Hit::new(0, (2, 2));

        assert!(Hit::new(0, (1, 1)).adjacent(&centre));
        assert!(Hit::new(0, (1, 2)).adjacent(&centre));
        assert!(Hit::new(0, (1, 3)).adjacent(&centre));
        assert!(Hit::new(0, (2, 1)).adjacent(&centre));
        assert!(Hit::new(0, (2, 3)).adjacent(&centre));
        assert!(Hit::new(0, (3, 1)).adjacent(&centre));
        assert!(Hit::new(0, (3, 2)).adjacent(&centre));
        assert!(Hit::new(0, (3, 3)).adjacent(&centre));

        assert!(Hit::new(0, (2, 2)).adjacent(&centre));

        assert!(!Hit::new(1, (0, 2)).adjacent(&centre));
        assert!(!Hit::new(0, (2, 0)).adjacent(&centre));
        assert!(!Hit::new(0, (2, 4)).adjacent(&centre));
        assert!(!Hit::new(0, (4, 2)).adjacent(&centre));
    }

    #[test]
    fn ccl_() {
        let hits = vec![
            Hit::new(0, (2, 2)),
            Hit::new(0, (3, 3)),
            Hit::new(0, (9, 9)),
            Hit::new(1, (3, 3)),
            Hit::new(1, (2, 2)),
            Hit::new(0, (4, 4)),
        ];

        assert_eq!(ccl(hits), vec![0, 0, 2, 3, 3, 0]);
    }
}
