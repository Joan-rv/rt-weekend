use glam::{UVec2, Vec3};
use image::{ImageReader, Rgb, RgbImage};
use indicatif::{ParallelProgressIterator, ProgressBar, ProgressStyle};
use minifb::{Window, WindowOptions};
use rayon::prelude::*;
use rt_weekend::rt::{
    BvhNode, Camera, CameraDesc, CheckerTexture, ColorTexture, Dielectric, Hittable, Lambertian,
    Metal, Sphere, render_pixel,
};
use std::env::args;
use std::sync::mpsc;
use std::thread;
use std::{f32::consts::PI, sync::Arc};

fn color_to_rgb(color: Vec3) -> Rgb<u8> {
    Rgb((256.0
        * color
            .map(|v| if v > 0.0 { v.sqrt() } else { v })
            .clamp(Vec3::splat(0.0), Vec3::splat(0.999)))
    .as_u8vec3()
    .to_array())
}

fn bouncing_balls_scene() -> (Camera, Box<dyn Hittable>) {
    let mut world: Vec<Arc<dyn Hittable>> = Vec::new();

    let checker = Arc::new(CheckerTexture {
        scale: 0.32,
        even: Arc::new(ColorTexture {
            albedo: Vec3::new(0.2, 0.3, 0.1),
        }),
        odd: Arc::new(ColorTexture {
            albedo: Vec3::new(0.9, 0.9, 0.9),
        }),
    });
    world.push(Arc::new(Sphere::stationary(
        Vec3::new(0.0, -1000.0, 0.0),
        1000.0,
        Arc::new(Lambertian { texture: checker }),
    )));

    for a in -11..11 {
        for b in -11..11 {
            let choose_mat: f32 = rand::random();
            let center = Vec3::new(
                a as f32 + 0.9 * rand::random::<f32>(),
                0.2,
                b as f32 + 0.9 * rand::random::<f32>(),
            );

            if (center - Vec3::new(4.0, 0.2, 0.0)).length() > 0.9 {
                if choose_mat < 0.8 {
                    let albedo = rand::random::<Vec3>() * rand::random::<Vec3>();
                    let center2 = center + Vec3::new(0.0, rand::random_range(0.0..0.5), 0.0);
                    world.push(Arc::new(Sphere::moving(
                        center,
                        center2,
                        0.2,
                        Arc::new(Lambertian::from_color(albedo)),
                    )));
                } else if choose_mat < 0.95 {
                    let albedo = rand::random::<Vec3>() * 0.5 + 0.5;
                    let fuzziness = rand::random_range(0.0..0.5);
                    world.push(Arc::new(Sphere::stationary(
                        center,
                        0.2,
                        Arc::new(Metal { albedo, fuzziness }),
                    )));
                } else {
                    world.push(Arc::new(Sphere::stationary(
                        center,
                        0.2,
                        Arc::new(Dielectric {
                            refraction_index: 1.5,
                        }),
                    )))
                }
            }
        }
    }

    world.push(Arc::new(Sphere::stationary(
        Vec3::new(0.0, 1.0, 0.0),
        1.0,
        Arc::new(Dielectric {
            refraction_index: 1.5,
        }),
    )));
    world.push(Arc::new(Sphere::stationary(
        Vec3::new(-4.0, 1.0, 0.0),
        1.0,
        Arc::new(Lambertian::from_color(Vec3::new(0.4, 0.2, 0.1))),
    )));
    world.push(Arc::new(Sphere::stationary(
        Vec3::new(4.0, 1.0, 0.0),
        1.0,
        Arc::new(Metal {
            albedo: Vec3::new(0.7, 0.6, 0.5),
            fuzziness: 0.0,
        }),
    )));

    let world = BvhNode::create(&mut world);

    let camera = Camera::new(&CameraDesc {
        aspect_ratio: 16.0 / 9.0,
        image_width: 400,
        samples_per_px: 500,
        max_depth: 50,

        fovy: 20.0 * PI / 180.0,
        look_from: Vec3::new(13.0, 2.0, 3.0),
        look_at: Vec3::new(0.0, 0.0, 0.0),

        defocus_angle: 0.6 * PI / 180.0,
        focus_dist: 10.0,
    });

    (camera, Box::new(world))
}

fn checkered_spheres_scene() -> (Camera, Box<dyn Hittable>) {
    let mut world: Vec<Arc<dyn Hittable>> = Vec::new();

    let checker = Arc::new(CheckerTexture {
        scale: 0.32,
        even: Arc::new(ColorTexture {
            albedo: Vec3::new(0.2, 0.3, 0.1),
        }),
        odd: Arc::new(ColorTexture {
            albedo: Vec3::new(0.9, 0.9, 0.9),
        }),
    });
    world.push(Arc::new(Sphere::stationary(
        Vec3::new(0.0, -10.0, 0.0),
        10.0,
        Arc::new(Lambertian {
            texture: checker.clone(),
        }),
    )));
    world.push(Arc::new(Sphere::stationary(
        Vec3::new(0.0, 10.0, 0.0),
        10.0,
        Arc::new(Lambertian { texture: checker }),
    )));

    let camera = Camera::new(&CameraDesc {
        aspect_ratio: 16.0 / 9.0,
        image_width: 1200,
        samples_per_px: 100,
        max_depth: 50,

        fovy: 20.0 * PI / 180.0,
        look_from: Vec3::new(13.0, 2.0, 3.0),
        look_at: Vec3::new(0.0, 0.0, 0.0),

        defocus_angle: 0.0 * PI / 180.0,
        focus_dist: 10.0,
    });

    let world = BvhNode::create(&mut world);

    (camera, Box::new(world))
}

fn earth_scene() -> anyhow::Result<(Camera, Box<dyn Hittable>)> {
    let earth_texture = Arc::new(
        ImageReader::open("res/earthmap10k.jpg")?
            .decode()?
            .into_rgb32f(),
    );
    let earth_surface = Arc::new(Lambertian {
        texture: earth_texture,
    });
    let globe = Box::new(Sphere::stationary(
        Vec3::new(0.0, 0.0, 0.0),
        2.0,
        earth_surface,
    ));

    let camera = Camera::new(&CameraDesc {
        aspect_ratio: 16.0 / 9.0,
        image_width: 1920,
        samples_per_px: 100,
        max_depth: 50,

        fovy: 20.0 * PI / 180.0,
        look_from: Vec3::new(0.0, 0.0, 12.0),
        look_at: Vec3::new(0.0, 0.0, 0.0),
        defocus_angle: 0.0,
        focus_dist: 2.0,
    });

    Ok((camera, globe))
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
    let scene_type = args().nth(1);
    let (camera, world) = match scene_type.as_ref().map(|s| s.as_str()) {
        Some(s) if s.starts_with("bouncing") => bouncing_balls_scene(),
        Some(s) if s.starts_with("earth") => earth_scene()?,
        _ => checkered_spheres_scene(),
    };

    let (sender, receiver) = mpsc::channel();

    let width = camera.width();
    let height = camera.height();
    let gfx_thread = thread::spawn(move || gfx_thread(width, height, receiver));

    let bar = ProgressBar::new((camera.width() * camera.height()).into());
    bar.set_style(
        ProgressStyle::with_template("{wide_bar} {pos}/{len} {elapsed_precise}/{duration_precise}")
            .unwrap(),
    );
    (0..camera.height())
        .into_par_iter()
        .flat_map(|y| (0..camera.width()).into_par_iter().map(move |x| (x, y)))
        .progress_with(bar)
        .try_for_each_with(sender, |sender, (x, y)| -> anyhow::Result<()> {
            let color = color_to_rgb(render_pixel(&camera, &*world, UVec2::new(x, y)));
            sender.send(Pixel { x, y, color })?;
            Ok(())
        })?;

    gfx_thread.join().unwrap()?;

    Ok(())
}
