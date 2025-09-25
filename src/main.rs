use glam::{UVec2, Vec3};
use image::{Rgb, RgbImage};
use std::{f32::consts::FRAC_PI_4, rc::Rc};

pub const INFINITY: f32 = std::f32::MAX;

pub struct Interval {
    pub min: f32,
    pub max: f32,
}

impl Interval {
    pub const EMPTY: Self = Self {
        min: INFINITY,
        max: -INFINITY,
    };
    pub const UNIVERSE: Self = Self {
        min: -INFINITY,
        max: INFINITY,
    };

    pub fn contains(&self, x: f32) -> bool {
        self.min <= x && x <= self.max
    }
    pub fn surrounds(&self, x: f32) -> bool {
        self.min < x && x < self.max
    }
}

pub fn random_unit_vec3() -> Vec3 {
    loop {
        let v = 2.0 * rand::random::<Vec3>() - 1.0;
        let len_sq = v.length_squared();
        if len_sq < 1.0
            && let Some(v) = v.try_normalize()
        {
            return v;
        }
    }
}

pub fn random_unit_disk_vec3() -> Vec3 {
    loop {
        let v = 2.0 * rand::random::<Vec3>() - 1.0;
        if v.length_squared() < 1.0 {
            return v;
        }
    }
}

pub struct Ray {
    pub origin: Vec3,
    pub direction: Vec3,
}

impl Ray {
    pub fn at(&self, t: f32) -> Vec3 {
        self.origin + self.direction * t
    }
}

pub struct HitRecord {
    pub point: Vec3,
    pub normal: Vec3,
    pub material: Rc<dyn Material>,
    pub t: f32,
    pub front_face: bool,
}

pub trait Hittable {
    fn hit(&self, ray: &Ray, valid_t: Interval) -> Option<HitRecord>;
}

pub trait Material {
    fn scatter(&self, ray: &Ray, record: &HitRecord) -> Option<(Vec3, Ray)>;
}

pub struct Camera {
    image_width: u32,
    image_height: u32,

    center: Vec3,
    pixel00: Vec3,
    pixel_du: Vec3,
    pixel_dv: Vec3,

    defocus_angle: f32,
    defocus_disk_u: Vec3,
    defocus_disk_v: Vec3,

    samples_per_px: u32,
    max_depth: u32,
}

pub struct CameraDesc {
    pub image_width: u32,
    pub aspect_raio: f32,

    pub look_from: Vec3,
    pub look_at: Vec3,

    pub defocus_angle: f32,
    pub focus_dist: f32,

    pub fovy: f32,

    pub samples_per_px: u32,
    pub max_depth: u32,
}

impl Camera {
    pub fn new(desc: &CameraDesc) -> Self {
        let image_height = (desc.image_width as f32 / desc.aspect_raio) as u32;
        let h = (desc.fovy / 2.0).tan();

        let viewport_height = 2.0 * h * desc.focus_dist;
        let viewport_width = viewport_height * desc.image_width as f32 / image_height as f32;

        let w = (desc.look_from - desc.look_at).normalize();
        let u = Vec3::Y.cross(w).normalize();
        let v = w.cross(u);

        let viewport_u = viewport_width * u;
        let viewport_v = viewport_height * -v;
        let pixel_du = viewport_u / (desc.image_width as f32);
        let pixel_dv = viewport_v / (image_height as f32);

        let viewport_top_left =
            desc.look_from - (desc.focus_dist * w) - viewport_u / 2.0 - viewport_v / 2.0;
        let pixel00 = viewport_top_left + 0.5 * (pixel_du + pixel_dv);

        let defocus_radius = desc.focus_dist * (desc.defocus_angle / 2.0).tan();
        let defocus_disk_u = u * defocus_radius;
        let defocus_disk_v = v * defocus_radius;
        Self {
            image_width: desc.image_width,
            image_height,
            center: desc.look_from,
            pixel00,
            pixel_du,
            pixel_dv,

            defocus_angle: desc.defocus_angle,
            defocus_disk_u,
            defocus_disk_v,

            samples_per_px: desc.samples_per_px,
            max_depth: desc.max_depth,
        }
    }

    pub fn width(&self) -> u32 {
        self.image_width
    }
    pub fn height(&self) -> u32 {
        self.image_height
    }

    fn get_ray(&self, coords: UVec2) -> Ray {
        let x = coords.x as f32 + rand::random::<f32>() - 0.5;
        let y = coords.y as f32 + rand::random::<f32>() - 0.5;
        let px_loc = self.pixel00 + x * self.pixel_du + y * self.pixel_dv;
        let origin = if self.defocus_angle > 0.0 {
            let offset = random_unit_disk_vec3();
            self.center + offset.x * self.defocus_disk_u + offset.y * self.defocus_disk_v
        } else {
            self.center
        };
        Ray {
            origin,
            direction: px_loc - origin,
        }
    }
}

fn ray_color(ray: &Ray, depth: u32, world: &impl Hittable) -> Vec3 {
    if depth == 0 {
        return Vec3::ZERO;
    }

    let pos_interval = Interval {
        min: 0.001,
        max: INFINITY,
    };
    if let Some(record) = world.hit(ray, pos_interval) {
        if let Some((attenuation, scattered)) = record.material.scatter(ray, &record) {
            return attenuation * ray_color(&scattered, depth - 1, world);
        }
        return Vec3::X;
    }

    let dir = ray.direction.normalize();
    let a = 0.5 * (dir.y + 1.0);
    (1.0 - a) * Vec3::new(1.0, 1.0, 1.0) + a * Vec3::new(0.5, 0.7, 1.0)
}

pub fn render_pixel(camera: &Camera, world: &impl Hittable, coords: UVec2) -> Vec3 {
    let mut color = Vec3::ZERO;
    for _ in 0..camera.samples_per_px {
        color += ray_color(&camera.get_ray(coords), camera.max_depth, world);
    }
    color / (camera.samples_per_px as f32)
}

struct Sphere {
    center: Vec3,
    radius: f32,
    material: Rc<dyn Material>,
}

impl Hittable for Sphere {
    fn hit(&self, ray: &Ray, valid_t: Interval) -> Option<HitRecord> {
        let oc = self.center - ray.origin;
        let a = ray.direction.length_squared();
        let h = ray.direction.dot(oc);
        let c = oc.length_squared() - self.radius * self.radius;
        let discriminant = h * h - a * c;
        if discriminant < 0.0 {
            return None;
        }

        let discriminant_sqrt = discriminant.sqrt();
        let mut root = (h - discriminant_sqrt) / a;
        if !valid_t.surrounds(root) {
            root = (h + discriminant_sqrt) / a;
            if !valid_t.surrounds(root) {
                return None;
            }
        }

        let point = ray.at(root);
        let mut normal = (point - self.center) / self.radius;
        let front_face = ray.direction.dot(normal) < 0.0;
        if !front_face {
            normal = -normal;
        }
        Some(HitRecord {
            point,
            normal,
            t: root,
            material: self.material.clone(),
            front_face,
        })
    }
}

pub struct Lambertian {
    pub albedo: Vec3,
}

impl Material for Lambertian {
    fn scatter(&self, _ray: &Ray, record: &HitRecord) -> Option<(Vec3, Ray)> {
        let mut out_dir = random_unit_vec3() + record.normal;
        if out_dir.length_squared() < 1e-4 {
            out_dir = record.normal;
        }
        Some((
            self.albedo,
            Ray {
                direction: out_dir,
                origin: record.point,
            },
        ))
    }
}

fn color_to_rgb(color: Vec3) -> Rgb<u8> {
    Rgb((256.0
        * color
            .map(|v| if v > 0.0 { v.sqrt() } else { v })
            .clamp(Vec3::splat(0.0), Vec3::splat(0.999)))
    .as_u8vec3()
    .to_array())
}

fn main() -> anyhow::Result<()> {
    let material = Rc::new(Lambertian {
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

    let mut image = RgbImage::new(camera.width(), camera.height());

    for y in 0..camera.height() {
        for x in 0..camera.width() {
            eprint!(
                "\r{}/{}",
                y * camera.width() + x,
                camera.width() * camera.height()
            );
            image.put_pixel(
                x,
                y,
                color_to_rgb(render_pixel(&camera, &world, UVec2::new(x, y))),
            );
        }
    }
    println!("");

    image.save("im.png")?;

    Ok(())
}
