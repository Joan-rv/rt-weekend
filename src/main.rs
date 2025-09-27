use glam::{UVec2, Vec3};
use image::{Rgb, RgbImage};
use indicatif::ProgressBar;
use minifb::{Window, WindowOptions};
use rayon::prelude::*;
use rt_weekend::rt::{
    Camera, CameraDesc, Dielectric, Hittable, Lambertian, Metal, Sphere, render_pixel,
};
use std::sync::mpsc;
use std::thread;
use std::{
    f32::consts::{FRAC_PI_4, PI},
    sync::Arc,
};

fn color_to_rgb(color: Vec3) -> Rgb<u8> {
    Rgb((256.0
        * color
            .map(|v| if v > 0.0 { v.sqrt() } else { v })
            .clamp(Vec3::splat(0.0), Vec3::splat(0.999)))
    .as_u8vec3()
    .to_array())
}

fn simple_scene() -> (Camera, Box<dyn Hittable>) {
    let material = Arc::new(Lambertian {
        albedo: Vec3::new(1.0, 0.0, 0.0),
    });
    let world = Sphere {
        center: Vec3::new(0.0, 0.0, 0.0),
        radius: 0.1,
        material,
    };
    let camera = Camera::new(&CameraDesc {
        image_width: 200,
        aspect_raio: 16.0 / 9.0,

        look_from: Vec3::new(0.0, 0.0, -1.0),
        look_at: Vec3::new(0.0, 0.0, 0.0),

        defocus_angle: 0.1,
        focus_dist: 1.0,
        fovy: FRAC_PI_4,

        max_depth: 32,
        samples_per_px: 10,
    });

    (camera, Box::new(world))
}

fn sample_scene() -> (Camera, Box<dyn Hittable>) {
    let camera = Camera::new(&CameraDesc {
        look_from: Vec3::new(-2.0, 2.0, 1.0),
        look_at: Vec3::new(0.0, 0.0, -1.0),
        defocus_angle: PI * 10.0 / 180.0,
        focus_dist: 3.4,
        image_width: 2000,
        aspect_raio: 16.0 / 9.0,
        fovy: PI * 20.0 / 180.0,
        samples_per_px: 100,
        max_depth: 50,
    });

    let mat_ground = Arc::new(Lambertian {
        albedo: Vec3::new(0.8, 0.8, 0.0),
    });
    let mat_center = Arc::new(Lambertian {
        albedo: Vec3::new(0.1, 0.2, 0.5),
    });
    let mat_left = Arc::new(Dielectric {
        refraction_index: 1.5,
    });
    let mat_bubble = Arc::new(Dielectric {
        refraction_index: 1.0 / 1.5,
    });
    let mat_right = Arc::new(Metal {
        albedo: Vec3::new(0.8, 0.6, 0.2),
        fuzziness: 1.0,
    });

    let spheres: Vec<_> = vec![
        Sphere {
            center: Vec3::new(0.0, 0.0, -1.0),
            radius: 0.5,
            material: mat_center.clone(),
        },
        Sphere {
            center: Vec3::new(0.0, -100.5, -1.0),
            radius: 100.0,
            material: mat_ground.clone(),
        },
        Sphere {
            center: Vec3::new(-1.0, 0.0, -1.0),
            radius: 0.5,
            material: mat_left.clone(),
        },
        Sphere {
            center: Vec3::new(-1.0, 0.0, -1.0),
            radius: 0.4,
            material: mat_bubble.clone(),
        },
        Sphere {
            center: Vec3::new(1.0, 0.0, -1.0),
            radius: 0.5,
            material: mat_right.clone(),
        },
    ]
    .into_iter()
    .map(|s| Box::new(s) as Box<dyn Hittable>)
    .collect();
    (camera, Box::new(spheres))
}

struct Pixel {
    x: u32,
    y: u32,
    color: Rgb<u8>,
}

fn gfx_thread(width: u32, height: u32, receiver: mpsc::Receiver<Pixel>) -> anyhow::Result<()> {
    let mut image = RgbImage::new(width, height);
    let mut window = Window::new(
        "rt",
        width as usize,
        height as usize,
        WindowOptions::default(),
    )?;
    window.set_target_fps(240);

    let mut buffer = vec![0x000000u32; width as usize * height as usize];
    'window: while window.is_open() {
        loop {
            match receiver.try_recv() {
                Ok(px) => {
                    image.put_pixel(px.x, px.y, px.color);
                    buffer[(px.y * width + px.x) as usize] =
                        u32::from_be_bytes([0x00, px.color[0], px.color[1], px.color[2]]);
                }
                Err(mpsc::TryRecvError::Disconnected) => break 'window,
                Err(mpsc::TryRecvError::Empty) => break,
            }
        }
        window.update_with_buffer(&buffer, width as usize, height as usize)?;
    }

    image.save("im.png")?;

    Ok(())
}

fn main() -> anyhow::Result<()> {
    let (camera, world) = sample_scene();

    let (sender, receiver) = mpsc::channel();

    let width = camera.width();
    let height = camera.height();
    let gfx_thread = thread::spawn(move || gfx_thread(width, height, receiver));

    let bar = ProgressBar::new((camera.width() * camera.height()).into());
    (0..camera.height())
        .into_par_iter()
        .flat_map(|y| (0..camera.width()).into_par_iter().map(move |x| (x, y)))
        .try_for_each_with(sender, |sender, (x, y)| -> anyhow::Result<()> {
            bar.inc(1);
            let color = color_to_rgb(render_pixel(&camera, &*world, UVec2::new(x, y)));
            sender.send(Pixel { x, y, color })?;
            Ok(())
        })?;
    bar.finish();

    gfx_thread.join().unwrap()?;

    Ok(())
}
