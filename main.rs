// Mandelbrot, multi thread.
#![allow(dead_code, unused_variables)]
extern crate crossbeam;


use std::str::FromStr;

/// Parse the string `s` as a coordinate pair, like `"400x600"` or `"1.0,0.5"`.
///
/// Specifically, `s` should have the form <left><sep><right>, where <sep> is
/// the character given by the `separator` argument, and <left> and <right> are both
/// strings that can be parsed by `T::from_str`.                                                         c
///
/// If `s` has the proper form, return `Some<(x, y)>`. If it doesn't parse
/// correctly, return `None`.
fn parse_pair<T: FromStr>(s: &str, separator: char) -> Option<(T, T)> {
    match s.find(separator) {
        None => None,
        Some(index) => match (T::from_str(&s[..index]), T::from_str(&s[index + 1..])) {
            (Ok(l), Ok(r)) => Some((l, r)),
            _ => None,
        },
    }
}

#[test]
fn test_parse_pair() {
    assert_eq!(parse_pair::<i32>("", ','), None);
    assert_eq!(parse_pair::<i32>("10,", ','), None);
    assert_eq!(parse_pair::<i32>(",10", ','), None);
    assert_eq!(parse_pair::<i32>("10,20", ','), Some((10, 20)));
    assert_eq!(parse_pair::<i32>("10,20xy", ','), None);
    assert_eq!(parse_pair::<f64>("0.5x", 'x'), None);
    assert_eq!(parse_pair::<f64>("0.5x1.5", 'x'), Some((0.5, 1.5)));
}

///	Parse	a	pair	of	floating-point	numbers	separated	by	a	comma	as	a	complex
///	number.
fn parse_complex(s: &str) -> Option<Complex<f64>> {
    match parse_pair(s, ',') {
        Some((re, im)) => Some(Complex { re, im }),
        None => None,
    }
}

#[test]
fn test_parse_complex() {
    assert_eq!(
        parse_complex("1.25,-0.0625"),
        Some(Complex {
            re: 1.25,
            im: -0.0625
        })
    );
    assert_eq!(parse_complex(",-0.0625"), None);
}

//  	The	following function	converts	from	image	space	to	complex	number	space:
/// Return the point on the complex plane corresponding to a given pixel in the
/// bitmap.
///
/// `bounds` is a pair giving the width and height of the bitmap. `pixel` is a
/// pair indicating a particular pixel in that bitmap. The `upper_left` and
/// `lower_right` parameters are points on the complex plane designating the
/// area our bitmap covers.
fn pixel_to_point(
    bounds: (usize, usize),
    pixel: (usize, usize),
    upper_left: (f64, f64),
    lower_right: (f64, f64),
) -> (f64, f64) {
    // It might be nicer to find the position of the *middle* of the pixel,
    // instead of its upper left corner, but this is easier to write tests for.
    let (width, height) = (lower_right.0 - upper_left.0, upper_left.1 - lower_right.1);
    (
        upper_left.0 + pixel.0 as f64 * width / bounds.0 as f64,
        upper_left.1 - pixel.1 as f64 * height / bounds.1 as f64,
    )
    // Why	subtraction	here?	pixel.1	increases	as	we	go	down,
    //	but	the	imaginary	component	increases	as	we	go	up.
}

#[test]
fn test_pixel_to_point() {
    assert_eq!(
        pixel_to_point((100, 100), (25, 75), (-1.0, 1.0), (1.0, -1.0)),
        (-0.5, -0.5)
    );
}

extern crate num;
use num::Complex;

/// Try to determine whether the complex number `c` is in the Mandelbrot set.
///
/// A number `c` is in the set if, starting with zero, repeatedly squaring and
/// adding `c` never causes the number to leave the circle of radius 2 centered
/// on the origin; the number instead orbits near the origin forever. (If the
/// number does leave the circle, it eventually flies away to infinity.)
///
/// If after `limit` iterations our number has still not left the circle, return
/// `None`; this is as close as we come to knowing that `c` is in the set.
///
/// If the number does leave the circle before we give up, return `Some(i)`, where
/// `i` is the number of iterations it took.
fn escape_time(c: Complex<f64>, limit: u32) -> Option<u32> {
    let mut z = Complex { re: 0.0, im: 0.0 };
    for i in 0..limit {
        z = z * z + c;
        if z.norm_sqr() > 4.0 {
            return Some(i);
        }
    }

    None
}

/// Render a rectangle of the Mandelbrot set into a buffer of pixels.
///
/// The `bounds` argument gives the width and height of the buffer `pixels`,
/// which holds one grayscale pixel per byte. The `upper_left` and `lower_right`
/// arguments specify points on the complex plane corresponding to the upper
/// left and lower right corners of the pixel buffer.
fn render(
    pixels: &mut [u8],
    bounds: (usize, usize),
    upper_left: (f64, f64),
    lower_right: (f64, f64),
) {
    assert!(pixels.len() == bounds.0 * bounds.1);

    for r in 0..bounds.1 {
        for c in 0..bounds.0 {
            let point = pixel_to_point(bounds, (c, r), upper_left, lower_right);
            pixels[r * bounds.0 + c] = match escape_time(
                Complex {
                    re: point.0,
                    im: point.1,
                },
                255,
            ) {
                None => 0,
                Some(count) => 255 - count as u8,
            };
            // If	escape_time	says	that	point	belongs	to	the	set,	render	colors	the corresponding	pixel	black	(0).
            // Otherwise,	render	assigns	darker	colors	to the	numbers	that	took	longer	to	escape	the	circle.
        }
    }
}

extern	crate	image;
use image::png::PNGEncoder;
use image::ColorType;
use std::fs::File;

// type	Result<T>	=	std::result::Result<T,	Error>  type alias

///	Write	the	buffer	`pixels`,	whose	dimensions	are	given	by	`bounds`,	to	the
///	file	named	`filename`.
fn write_image(
    filename: &str,
    pixels: &[u8],
    bounds: (usize, usize),
) -> Result<(), std::io::Error> {
    let output = File::create(filename)?;

    let encoder = PNGEncoder::new(output);
    encoder.encode(
        &pixels,
        bounds.0 as u32,
        bounds.1 as u32,
        ColorType::Gray(8),
    )?;

    Ok(())
}

use std::io::Write;

fn main() {
    let args: Vec<String> = std::env::args().collect();

    if args.len() != 5 {
        writeln!(
            std::io::stderr(),
            "Usage: mandelbrot FILE PIXELS UPPERLEFT LOWERRIGHT"
        )
        .unwrap();
        writeln!(
            std::io::stderr(),
            "Example: {} mandel.png 1000x750 -1.20,0.35 -1,0.20",
            args[0]
        )
        .unwrap();
        std::process::exit(1);
    }

    let bounds = parse_pair(&args[2], 'x').expect("error parsing image dimensions");
    let upper_left = parse_pair(&args[3], ',').expect("error parsing upper left corner point");
    let lower_right = parse_pair(&args[4], ',').expect("error parsing lower right corner point");

    // After collecting	the	command-line	arguments	into	a	vector	of	Strings
    // we parse	each	one	and	then	begin	calculations.

    let mut pixels = vec![0; bounds.0 * bounds.1];
    // A	macro	call	vec![v;	n]	creates	a	vector	n	elements	long	whose	elements are	initialized	to	v,
    // so	the	preceding	code	creates	a	vector	of	zeros	whose length	is	bounds.0	*	bounds.1,
    // where	bounds	is	the	image	resolution parsed	from	the	command	line.
    // We’ll	use	this	vector	as	a	rectangular array	of	one-byte	grayscale	pixel	values

    // Then this	calls	the	render	function	to	actually	compute	the	image
    //render(&mut pixels[..], bounds, upper_left, lower_right);
    // the above single line calling render() is replaced in the multi threaded version by the code below
    
    let threads = 8;                               // here we decide to use eight threads.
    let rows_per_band = bounds.1 / threads + 1;    // how many rows of pixels each band should have
    
    {
		let bands: Vec<&mut[u8]> =
		    pixels.chunks_mut(rows_per_band*bounds.0).collect();     //  divide pixel buffer into bands
		    
		crossbeam::scope( | spawner | {   // crossbeam::scope calls the closure passing as the spawner
                                          //argument a value the closure can use to create new threads.
                                          
                                   // The crossbeam::scope function waits for all such threads to finish execution before returning itself
			
			for (i, band) in bands.into_iter().enumerate() {  // Here we iterate over the pixel buffer’s bands. The into_iter() iterator gives
                                                              // each iteration of the loop body exclusive ownership of one band,
                                                              // ensuring that only one thread can write to it at a time.
                let top = rows_per_band * i;
                let height = band.len() / bounds.0;
                let band_bounds = (bounds.0, height);
                let band_upper_left =
                     pixel_to_point(bounds, (0, top), upper_left, lower_right);
                let band_lower_right =
                     pixel_to_point(bounds, (bounds.0, top + height),
                                       upper_left, lower_right);

                spawner.spawn(move |_| {                      // we create a thread, running the closure move || { ... }.
                      render(band, band_bounds, band_upper_left, band_lower_right); // only the closure may use the mutable slice 'band'.
                });
            }
        });
    }
 
    // crossbeam::scope call ensures that all threads have completed before it returns, meaning that it is safe to save
    // the image to a file, which is our next action. 
    
    write_image(&args[1], &pixels, bounds).expect("error	writing	PNG	file");
}










