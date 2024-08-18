use std::fs;
use std::io::Write;

const XCONST: (f64, f64, f64) = (0.0, -90.0, 0.0);
const YCONST: (f64, f64, f64) = (-90.0, -90.0, 0.0);
const ZCONST: (f64, f64, f64) = (0.0, 0.0, 0.0);

pub fn create_box_detector(
    layers: &Vec<f64>,
    filename: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut id: u64 = 0;
    let mut file = fs::File::create(filename)?;

    let rotations = [XCONST, XCONST, YCONST, YCONST, ZCONST, ZCONST];
    for radius in layers {
        let dev = radius / 2.0;
        let translations = [
            (-dev, -dev, -dev),
            (dev, -dev, -dev),
            (-dev, -dev, -dev),
            (-dev, dev, -dev),
            (-dev, -dev, -dev),
            (-dev, -dev, dev),
        ];
        let radius_u = radius.round() as u64;
        for (transl, rot) in translations.iter().zip(rotations) {
            file.write_fmt(format_args!(
                "{} {} {} {} {} {} {} {} {} {} {}\n",
                id,
                radius,
                radius,
                transl.0,
                transl.1,
                transl.2,
                radius_u,
                radius_u,
                rot.0,
                rot.1,
                rot.2
            ))?;
            id += 1;
        }
    }
    Ok(())
}
