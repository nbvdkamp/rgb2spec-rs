use rgb2spec::optimize::{gamut::Gamut, optimize};

fn main() {
    let args: Vec<String> = std::env::args().collect();

    if args.len() < 3 {
        println!(
            "Syntax: rgb2spec_opt <resolution> <output> [<gamut>]
                  where <gamut> is one of sRGB,eRGB,XYZ,ProPhotoRGB,ACES2065_1,REC2020"
        );
        std::process::exit(-1);
    }

    let gamut = if args.len() > 3 {
        match Gamut::parse(&args[3]) {
            Some(gamut) => gamut,
            None => {
                eprintln!("Could not parse gamut `{}'!", args[3]);
                std::process::exit(-1);
            }
        }
    } else {
        Gamut::SRGB
    };

    let res: usize = match args[1].parse() {
        Ok(res) => res,
        Err(_) => {
            eprintln!("Could not parse resolution `{}'!", args[1]);
            std::process::exit(-1);
        }
    };

    if res > 0xFFFF {
        eprintln!("Resolution too large `{}'!", args[1]);
        std::process::exit(-1);
    }

    print!("Optimizing spectra ...");
    let model = match optimize(gamut, res) {
        Ok(m) => m,
        Err(e) => {
            let hint = if res < 10 {
                "\nTry a resolution of at least 10."
            } else {
                ""
            };

            eprintln!("\n An error occured during optimizing:\n{e}{hint}");
            std::process::exit(-1);
        }
    };
    println!(" done.");

    match model.save(&args[2]) {
        Ok(()) => (),
        Err(e) => eprintln!("Saving file failed: {e}"),
    }
}
